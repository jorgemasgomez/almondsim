# Prepare only: running this file explicitly starts the speed breeding scenario.
.libPaths(c(Sys.getenv('ALMONDSIM_R_LIBRARY',
  'C:/Users/franc/Documents/Codex/2026-10-04/es/work/r-library'),
  'C:/Users/franc/AppData/Local/R/win-library/4.4',.libPaths()))
files<-Filter(Negate(is.null),lapply(sys.frames(),function(frame) frame$ofile))
flag<-grep('^--file=',commandArgs(),value=TRUE)
script<-if(length(files)) tail(files,1)[[1]] else sub('^--file=','',flag[1])
folder<-dirname(normalizePath(script,winslash='/',mustWork=TRUE));root<-dirname(folder)
settings_env<-new.env(parent=baseenv())
sys.source(file.path(folder,'settings.R'),envir=settings_env)
settings<-settings_env$settings
source(file.path(root,'redesigned_schemes.R'))
source(file.path(root,'00_Burn_in','GlobalParameters.R'))
source(file.path(root,'compatible_crosses.R'))
source(file.path(root,'ocs_discrete.R'))
source(file.path(root,'speed_breeding.R'))
source(file.path(root,'optisel_adapter.R'))
base_run<-settings$base_run
configuration<-readRDS(file.path(base_run,'configuration.rds'))
run<-tempfile(paste0(settings$id,'_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),
  tmpdir=file.path(redesign_data_dir(root,configuration$config),'redesigned_runs'))
dir.create(run,recursive=TRUE)
saveRDS(list(settings=settings,source_run=base_run),file.path(run,'configuration.rds'))
RNGkind("L'Ecuyer-CMRG");set.seed(configuration$config$seed);stream<-.Random.seed
started<-proc.time()[['elapsed']];outputs<-list()
for(rep in seq_len(settings$reps)) {
  directory<-file.path(run,sprintf('rep_%03d',rep));dir.create(directory)
  load(file.path(base_run,sprintf('rep_%03d',rep),'burnin_redesigned.RData'),envir=.GlobalEnv)
  SP$nThreads<-6L
  state$config$parent_stage[['2']]<-'Juvenile2'
  state$config$speed_breeding<-TRUE
  state$config$early_selection<-TRUE;state$config$progeny<-75L
  state$config$seedlings_keep<-500L;state$config$entrants<-8L
  state$config$ocs_alpha<-settings$ocs_alpha
  state$config$opr_max_diversity_loss<-NA_real_
  state$config$scenario_id<-settings$id
  # Reinstall from clean functions for each paired replicate.
  source(file.path(root,'redesigned_schemes.R'))
  source(file.path(root,'compatible_crosses.R'))
  source(file.path(root,'ocs_discrete.R'))
  speed_install(file.path(root,'optisel_adapter.R'))
  assign('.Random.seed',parallel::nextRNGStream(stream),.GlobalEnv)
  records<-list();log<-file.path(directory,'execution_speed.log')
  for(k in seq_len(40L)) {
    state$year<-40L+k
    jsonlite::write_json(list(status='running',rep=rep,scenario=settings$id,year=state$year),
      file.path(directory,'progress_speed.json'),auto_unbox=TRUE,pretty=TRUE)
    begin<-proc.time()[['elapsed']]
    row<-tryCatch(redesign_year(state,'GS',2L,rep),error=function(e) {
      writeLines(conditionMessage(e),file.path(directory,'error_speed.txt'))
      jsonlite::write_json(list(status='failed',rep=rep,year=state$year,error=conditionMessage(e)),
        file.path(directory,'progress_speed.json'),auto_unbox=TRUE,pretty=TRUE)
      stop(e)
    })
    records[[k]]<-row
    cat(sprintf('Year %d completed in %.2f seconds\n',state$year,proc.time()[['elapsed']]-begin),file=log,append=TRUE)
    saveRDS(do.call(rbind,records),file.path(directory,'results_partial.rds'))
    saveRDS(list(crosses=state$speed_cross_trace,ocs=tail(state$ocs_trace,1),solver=state$speed_ocs_solver),
      file.path(directory,sprintf('cross_trace_year_%03d.rds',state$year)))
  }
  result<-do.call(rbind,records)
  stopifnot(nrow(result)==40L,all(result$meanMotherAge>=4),all(result$meanFatherAge>=2))
  saveRDS(result,file.path(directory,paste0('results_',settings$id,'.rds')))
  jsonlite::write_json(list(status='complete',rep=rep,year=80L),
    file.path(directory,'progress_speed.json'),auto_unbox=TRUE,pretty=TRUE)
  outputs[[rep]]<-result;stream<-parallel::nextRNGStream(stream)
}
combined<-do.call(rbind,outputs)
stopifnot(nrow(combined)==40L*settings$reps,!anyDuplicated(combined[c('rep','scenario','year')]))
saveRDS(combined,file.path(run,'results_all_schemes.rds'))
saveRDS(list(complete=TRUE,replicates=settings$reps),file.path(run,'status.rds'))
saveRDS(list(elapsed_seconds=proc.time()[['elapsed']]-started,completed=Sys.time()),file.path(run,'total_runtime.rds'))
cat('COMPLETE: ',run,'\n',sep='')
