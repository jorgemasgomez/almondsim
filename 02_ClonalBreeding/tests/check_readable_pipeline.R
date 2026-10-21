# Regression checks for the source reorganization. No founder bank, ASReml
# licence or breeding trajectories are required; no production data are edited.
# Run from the repository root: Rscript 02_ClonalBreeding/tests/check_readable_pipeline.R
root <- normalizePath('02_ClonalBreeding',winslash='/',mustWork=TRUE)
engine <- new.env(parent=.GlobalEnv)
sys.source(file.path(root,'redesigned_schemes.R'),envir=engine)

# Compare executable syntax, not whitespace/comments or source-file locations.
# Scientific signatures were captured before refactoring. The controller's
# expected signature permits exactly the replacement of its old copy block by
# freeze_pipeline_source(original_root,root); all its other instructions match.
fingerprint <- function(fun) {
  path <- tempfile();on.exit(unlink(path))
  syntax <- parse(text=deparse(body(fun),width.cutoff=500L),keep.source=FALSE)[[1L]]
  saveRDS(list(formals=formals(fun),body=syntax),path,version=2,compress=FALSE)
  unname(tools::md5sum(path))
}
reference <- read.csv(file.path(root,'tests','pipeline_reference_signatures.csv'))
actual <- vapply(reference$function_name,function(name)fingerprint(engine[[name]]),character(1))
stopifnot(identical(unname(actual),reference$ast_md5))
cat('PASS:',nrow(reference)-1L,'unchanged scientific/runtime signatures and the controller with only its source-freeze update.\n')

# There must be a single explicit definition per module-level function. PT-only
# and repeated-record variants have distinct names instead of late overrides.
modules <- list.files(file.path(root,'R'),pattern='[.]R$',full.names=TRUE)
definitions <- unlist(lapply(modules,function(file) {
  expressions <- parse(file)
  unlist(lapply(expressions,function(e) {
    if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.call(e[[3]])&&
       identical(e[[3]][[1]],as.name('function')))as.character(e[[2]]) else NULL
  }))
}))
stopifnot(!anyDuplicated(definitions))
catalogue <- engine$redesign_all_scenarios(root)
aliases <- engine$scheme_aliases()
stopifnot(nrow(catalogue)==18L,length(aliases)==18L,
          setequal(unname(aliases),catalogue$scenario_id))
for(alias in names(aliases))stopifnot(identical(
  engine$select_breeding_schemes(root,alias)$scenario_id,unname(aliases[[alias]])))
stopifnot(inherits(try(engine$select_breeding_schemes(root,'not-a-scheme'),silent=TRUE),'try-error'))
cat('PASS: unique definitions, 18 catalogue entries, and all manuscript aliases.\n')

# Validation/help must not execute an experiment, and valid CLI options must
# reach the existing controller with the requested schemes and integer settings.
spy <- new.env(parent=engine)
spy$run_breeding_experiment <- function(root,schemes,config)list(root=root,schemes=schemes,config=config)
cli <- engine$pipeline_cli;environment(cli)<-spy
invisible(capture.output(stopifnot(is.null(cli(root,args='--help')))))
request <- cli(root,single=TRUE,args=c('--scheme=GS6-E','--replicates=2','--workers=1','--threads=3'))
stopifnot(identical(request$schemes,'GS6-E'),request$config$reps==2L,
          request$config$workers==1L,request$config$threads==3L)
for(args in list('--scheme=all',c('--scheme=GS6','--replicates=0'),
                c('--scheme=GS6','--unknown=1'),c('--scheme=GS6','--scheme=GS9')))
  stopifnot(inherits(try(cli(root,single=TRUE,args=args),silent=TRUE),'try-error'))
cat('PASS: CLI validation/dispatch without launching simulations.\n')

# Nested sourcing must resolve the launcher itself, not the caller's filename.
# Checking an isolated environment also catches accidental global-only loading.
for (launcher in c('run_scheme.R','run_experiment.R')) {
  runner <- new.env(parent=.GlobalEnv)
  sys.source(file.path(root,launcher),envir=runner)
  stopifnot(exists('run_breeding_scheme',envir=runner,inherits=FALSE),
            exists('run_breeding_experiment',envir=runner,inherits=FALSE),
            identical(runner$.pipeline_root,root))
}
cat('PASS: both launchers can be sourced from another script without starting an experiment.\n')

# A worker must load the copied modules without relying on the original checkout.
frozen <- tempfile('almondsim_frozen_');dir.create(frozen)
engine$freeze_pipeline_source(root,frozen)
snapshot <- new.env(parent=.GlobalEnv)
sys.source(file.path(frozen,'redesigned_schemes.R'),envir=snapshot)
stopifnot(identical(snapshot$redesign_config(),engine$redesign_config()))
manifest <- read.csv(file.path(frozen,'source_manifest.csv'))
stopifnot(identical(unname(tools::md5sum(file.path(frozen,manifest$file))),manifest$md5))
worker_script <- file.path(frozen,'check_frozen_worker.R')
writeLines(c(
  paste0('.libPaths(',paste(capture.output(dput(.libPaths())),collapse=''),')'),
  paste0('directory <- ',encodeString(frozen,quote='"')),
  "source(file.path(directory,'redesigned_schemes.R'))",
  "stopifnot(nrow(redesign_all_scenarios())==18L)",
  paste0("stopifnot(length(list.files(file.path(directory,'R')))==",length(modules),"L)"),
  "cat('Frozen worker module load passed.\\n')"),worker_script)
# A clean Rscript process exercises the same frozen-code loading boundary
# without requiring a local socket or starting a simulation worker job.
rscript <- file.path(R.home('bin'),'Rscript.exe')
if(!file.exists(rscript))rscript <- file.path(R.home('bin'),'Rscript')
output <- system2(rscript,shQuote(worker_script),stdout=TRUE,stderr=TRUE)
status <- attr(output,'status')
stopifnot(is.null(status)||status==0L,any(grepl('Frozen worker module load passed.',output,fixed=TRUE)))
unlink(frozen,recursive=TRUE)
cat('PASS: frozen source checksums and clean-process worker loading.\n')
