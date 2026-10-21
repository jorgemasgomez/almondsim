# Provide timing, storage checks, diagnostics and failure isolation. Timing is inclusive: never sum a parent task with its nested tasks.
# Public function names are retained for existing launchers and saved workflows.

# Capture an individual failure as an outcome so unrelated scenarios can continue.
redesign_attempt_scenario <- function(action) {
 tryCatch(list(ok=TRUE,result=action()),error=function(e)
   list(ok=FALSE,error=conditionMessage(e)))
}

# Resolve external output storage and refuse paths inside the source repository.
redesign_data_dir <- function(root,config=list()) {
  configured <- config$data_dir
  if(is.null(configured) || !nzchar(configured)) configured <- Sys.getenv('ALMONDSIM_DATA_DIR',unset='')
  if(!nzchar(configured)) configured <- file.path(dirname(dirname(normalizePath(root,winslash='/',mustWork=TRUE))),'almondsim_data')
  if(!grepl('^([A-Za-z]:[/\\\\]|/)',configured)) configured <- file.path(getwd(),configured)
  resolved <- normalizePath(configured,winslash='/',mustWork=FALSE)
  source_root <- normalizePath(root,winslash='/',mustWork=TRUE)
  repository <- if(basename(source_root)=='02_ClonalBreeding') dirname(source_root) else source_root
  if(tolower(resolved)==tolower(repository)||startsWith(tolower(resolved),paste0(tolower(repository),'/'))) {
    stop('Simulation data must be stored outside the source repository.')
  }
  resolved
}

# Wrap a task with inclusive timing and diagnostics, preserving the returned result or raised error.
redesign_timer <- function(original,label) {
  force(original);force(label)
  function(...) {
    started <- Sys.time(); cpu <- proc.time(); success <- FALSE
    on.exit({
      destination <- getOption('almond.timing_file')
      if(!is.null(destination)) {
        duration <- proc.time()-cpu
        row <- data.frame(rep=getOption('almond.timing_rep'),
          scenario=getOption('almond.timing_scenario'),year=getOption('almond.timing_year'),
          task=label,started=format(started,'%Y-%m-%dT%H:%M:%S%z'),
          elapsed_seconds=unname(duration['elapsed']),
          cpu_seconds=unname(duration['user.self']+duration['sys.self']),success=success)
        memory <- redesign_memory()
        row$working_set_gb <- memory$working_set_gb
        row$private_gb <- memory$private_gb
        row$peak_private_gb <- memory$peak_private_gb
        tryCatch(write.table(row,destination,sep='\t',row.names=FALSE,
          col.names=!file.exists(destination),append=file.exists(destination),quote=FALSE),
          error=function(e) cat('Timing log could not be written:',conditionMessage(e),'\n'))
      }
    },add=TRUE)
    result <- withCallingHandlers(original(...),
      warning=function(w) redesign_diagnostic(label,'warning',w),
      error=function(e) redesign_diagnostic(label,'error',e))
    success <- TRUE;result
  }
}

# Check the licensed pedigree solver in a separate process before scheduling pedigree fits.
redesign_check_asreml <- function() {
  # A separate preflight process releases its license before pedigree starts.
  script <- tempfile(fileext='.R');on.exit(unlink(script),add=TRUE)
  libraries <- paste(capture.output(dput(.libPaths())),collapse='')
  writeLines(c(paste0('.libPaths(',libraries,')'),
    "library(asreml)",
    "d <- data.frame(y=c(1,2,3,5,4,6),g=factor(rep(1:3,2)))",
    "invisible(asreml::asreml(y~1,random=~g,data=d,trace=FALSE))"),script)
  code <- system2(file.path(R.home('bin'),'Rscript.exe'),shQuote(script))
  if(code!=0L) stop('ASReml preflight failed before starting the experiment.')
  invisible(TRUE)
}

# Read process memory metrics when available; return missing metrics if the optional ps package is absent.
redesign_memory <- function() {
  result <- list(working_set_gb=NA_real_,private_gb=NA_real_,peak_private_gb=NA_real_)
  if(requireNamespace('ps',quietly=TRUE)) {
    memory <- tryCatch(ps::ps_memory_info(),error=function(e) NULL)
    if(!is.null(memory)) {
      result$working_set_gb <- unname(memory['rss'])/1024^3
      result$private_gb <- unname(memory['mem_private'])/1024^3
      result$peak_private_gb <- unname(memory['peak_pagefile'])/1024^3
    }
  }
  result
}

# Write a condition, call stack and memory snapshot to the configured diagnostic log.
redesign_diagnostic <- function(task,kind,condition) {
  destination <- getOption('almond.diagnostic_file')
  if(is.null(destination)) return(invisible(NULL))
  report <- c(paste(format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),kind,task),
    conditionMessage(condition),paste('Call:',paste(deparse(conditionCall(condition)),collapse=' ')),
    paste('Memory:',paste(names(redesign_memory()),unlist(redesign_memory()),collapse='; ')),
    'Call stack:',vapply(sys.calls(),function(x) paste(deparse(x),collapse=' '),character(1)))
  tryCatch(cat(paste(report,collapse='\n'),'\n',file=destination,append=TRUE),
    error=function(e) cat('Diagnostic:',conditionMessage(condition),'\n'))
  invisible(NULL)
}

