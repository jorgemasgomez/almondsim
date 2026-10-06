# Experimental scheme only; the main comparison remains at 25 schemes.
.libPaths(c('C:/Users/franc/AppData/Local/R/win-library/4.4',.libPaths()))
script_files <- Filter(Negate(is.null),lapply(sys.frames(),function(frame) frame$ofile))
flag <- grep('^--file=',commandArgs(),value=TRUE)
script <- if(length(script_files)) tail(script_files,1)[[1]] else sub('^--file=','',flag[1])
folder <- dirname(normalizePath(script,winslash='/',mustWork=TRUE))
root <- dirname(folder)
base_run <- 'C:/Users/franc/Desktop/Proyectos/Almond_sim/almondsim_data/redesigned_runs/run_20261005_202329_6f1453a95184'
source(file.path(root,'redesigned_schemes.R'))
cfg <- readRDS(file.path(base_run,'configuration.rds'))
config <- cfg$config;config$reps <- 4L;config$threads <- 6L
base <- cfg$scenarios[cfg$scenarios$scenario_id=='GS_Y2_Early1500_OPR_Gain',,drop=FALSE]
rare <- base
rare$scenario_id <- 'GS_Y2_Early1500_RareAlleleExchange'
rare$folder_name <- '26_GS_2_Early1500_RareAlleleExchange'
rare$folder <- folder
run <- tempfile(paste0('rare_exchange_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),
  tmpdir=file.path(redesign_data_dir(root,config),'redesigned_runs'))
dir.create(run,recursive=TRUE)
saveRDS(list(config=config,scenarios=rbind(base,rare),source_run=base_run,
  rare_markers=1000L,gain_weight=.5,protect_last_copy=TRUE,
  max_new_parents=8L,marker_importance_weight=.5,marker_effect_weight=.5),file.path(run,'configuration.rds'))
RNGkind("L'Ecuyer-CMRG");set.seed(config$seed);stream <- .Random.seed
started <- proc.time()[['elapsed']]
for(rep in 1:4) {
  source(file.path(root,'redesigned_schemes.R'),local=.GlobalEnv)
  worker <- redesign_worker
  code <- deparse(body(worker),width.cutoff=500L)
  index <- grep('source\\(file.path\\(job\\$root, "opr_genomic.R"\\)',code)
  stopifnot(length(index)==1L)
  code <- append(code,c(sprintf('source(%s,local=.GlobalEnv)',deparse(file.path(folder,'rare_allele_replacement.R')))),after=index)
  body(worker) <- parse(text=paste(code,collapse='\n'))[[1]]
  dir <- file.path(run,sprintf('rep_%03d',rep));dir.create(dir)
  stopifnot(file.copy(file.path(base_run,sprintf('rep_%03d',rep),'burnin_redesigned.RData'),dir))
  options(rare.output_dir=dir)
  job <- list(rep=rep,root=root,config=config,scenarios=rare,methods='GS',phase='experimental',
    seed=stream,libraries=.libPaths(),directory=dir)
  result <- worker(job)
  if(!isTRUE(result$ok)) stop('Experimental scheme failed: ',result$error)
  result$result$oprMaxDiversityLoss <- NA_real_
  result$result$rareMarkerCount <- 1000L
  result$result$rareGainWeight <- .5
  saveRDS(result$result,file.path(dir,paste0('results_',rare$scenario_id,'.rds')))
  baseline_file <- file.path(base_run,sprintf('rep_%03d',rep),paste0('results_',base$scenario_id,'.rds'))
  stopifnot(file.copy(baseline_file,dir))
  cat('Completed experimental replicate ',rep,'/4\n',sep='')
  stream <- parallel::nextRNGStream(stream)
}
frames <- unlist(lapply(1:4,function(rep) lapply(
  file.path(run,sprintf('rep_%03d',rep),paste0('results_',c(base$scenario_id,rare$scenario_id),'.rds')),readRDS)),recursive=FALSE)
columns <- unique(unlist(lapply(frames,names)))
frames <- lapply(frames,function(x) { for(n in setdiff(columns,names(x))) x[[n]] <- NA; x[columns] })
combined <- do.call(rbind,frames)
stopifnot(nrow(combined)==320L,!anyDuplicated(combined[c('rep','scenario','year')]))
saveRDS(combined,file.path(run,'results_all_schemes.rds'))
saveRDS(list(complete=TRUE,replicates=4L,scenarios=rbind(base,rare)),file.path(run,'status.rds'))
saveRDS(list(elapsed_seconds=proc.time()[['elapsed']]-started,completed=Sys.time()),file.path(run,'total_runtime.rds'))
cat('COMPLETE: ',run,'\n',sep='')
