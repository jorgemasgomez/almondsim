# Fill and burn in one replicate, then run its selected scenarios from paired starting states. A failed scenario is recorded; other scenarios continue.
# Public function names are retained for existing launchers and saved workflows.

# Execute one replicate job, retaining partial results, error messages and frozen-code provenance.
redesign_worker <- function(job) {
 .libPaths(job$libraries)
 options(almond.input_dir=job$root,almond.nThreads=job$config$threads)
 source(file.path(job$root,'redesigned_schemes.R'),local=.GlobalEnv)
 assign('.Random.seed',job$seed,envir=.GlobalEnv)
 job$config$environment_seed <- job$config$seed+900000L+as.integer(job$rep)
 dir.create(job$directory,recursive=TRUE,showWarnings=FALSE)
 timing_file<-file.path(job$directory,paste0('timings_',job$phase,'.tsv'))
 options(almond.timing_file=timing_file,almond.timing_rep=job$rep,
   almond.timing_scenario='initialization',almond.timing_year=0L,
   almond.diagnostic_file=file.path(job$directory,paste0('diagnostics_',job$phase,'.log')))
 install_timers<-function() {
   for(name in c('redesign_setup','redesign_advance','redesign_measure','redesign_fit',
     'redesign_predict','redesign_choose','redesign_renew','redesign_ocs',
     'randCrossGamSI','redesign_record','redesign_year','redesign_allocate','redesign_measure_parents')) {
     if(exists(name,envir=.GlobalEnv,inherits=FALSE))
       assign(name,redesign_timer(get(name,envir=.GlobalEnv),name),envir=.GlobalEnv)
   }
 }
 status_file<-file.path(job$directory,paste0('progress_',job$phase,'.json'))
 progress<-function(status,scenario,year=0L,error=NULL) {
   options(almond.timing_scenario=scenario,almond.timing_year=year)
   jsonlite::write_json(list(status=status,rep=job$rep,phase=job$phase,
     scenario=scenario,year=year,memory=redesign_memory(),
     updated=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),error=error),
     status_file,pretty=TRUE,auto_unbox=TRUE)
 }
 log<-file(file.path(job$directory,paste0('execution_',job$phase,'.log')),'wt')
 sink(log);on.exit({sink();close(log)},add=TRUE)
 base_file<-file.path(job$directory,'burnin_redesigned.RData')
 initialized<-redesign_attempt_scenario(function() {
   progress('running','initialization')
   if(job$phase=='parallel') {
     install_timers()
     state<-redesign_setup(job$root,job$config)
     burnin<-vector('list',job$config$burnin)
     for(year in seq_len(job$config$burnin)) {
       state$year<-year;progress('running','Burn_in',year)
       cat('Burn-in year',year,'\n')
       burnin[[year]]<-redesign_year(state,'Pheno',9L,job$rep)
       burnin[[year]]$scenario<-'Burn_in'
     }
     save(state,burnin,SP,locuscompt_position_vector,self_compatible_allele,file=base_file)
     saveRDS(redesign_bind(burnin),file.path(job$directory,'burnin_results.rds'))
     saveRDS(state$environment_p,file.path(job$directory,'annual_environment.rds'))
   } else if(!file.exists(base_file)) stop('This replicate has no successfully completed burn-in.')
   TRUE
 })
 if(!initialized$ok) {
   progress('failed','Burn_in',getOption('almond.timing_year'),initialized$error)
   writeLines(initialized$error,file.path(job$directory,paste0('error_',job$phase,'.txt')))
   status<-list(ok=FALSE,rep=job$rep,error=initialized$error,fatal_initialization=TRUE)
   saveRDS(status,file.path(job$directory,paste0('worker_status_',job$phase,'.rds')))
   return(status)
 }
 outputs<-list();statuses<-list()
 selected<-job$scenarios[job$scenarios$method%in%job$methods,,drop=FALSE]
 for(index in seq_len(nrow(selected))) {
   scenario<-selected$scenario_id[index]
   records<-list();completed_year<-job$config$burnin
   started<-Sys.time()
   outcome<-redesign_attempt_scenario(function() {
     # Reinstall clean functions and load a fresh paired starting population.
     # A failed scenario must not contaminate the next scenario's state or RNG.
     source(file.path(job$root,'redesigned_schemes.R'),local=.GlobalEnv)
     source(file.path(job$root,'ocs_discrete.R'),local=.GlobalEnv)
     source(file.path(job$root,'optimize_contributions.R'),local=.GlobalEnv)
     source(file.path(job$root,'compatible_crosses.R'),local=.GlobalEnv)
     load(base_file,envir=.GlobalEnv)
     state<-get('state',.GlobalEnv)
     if(!identical(state$config$design_version,job$config$design_version)) stop('Burn-in source design differs from the current implementation.')
     if(is.null(state$stage_history)) stop('Burn-in lacks stage-specific training archives.')
     SP$nThreads<-job$config$threads
     assign('.Random.seed',parallel::nextRNGStream(job$seed),envir=.GlobalEnv)
     method<-selected$method[index];years<-selected$first_parental_use[index]
     state$config$ocs_enabled<-selected$ocs_enabled[index]
     state$config$ocs_target_Ne<-selected$ocs_target_Ne[index]
     state$config$early_selection<-selected$early_selection[index]
     state$config$progeny<-selected$progeny_per_cross[index]
     state$config$seedlings_keep<-selected$seedlings_keep[index]
     state$config$entrants<-if(is.na(selected$renewal_percent[index]))
       selected$new_parents_per_year[index] else as.integer(state$config$parents*selected$renewal_percent[index]/100)
     state$config$scenario_id<-scenario
     state$config$speed_breeding<-isTRUE(selected$speed_breeding[index])
     state$config$minimum_pollen_age<-if(state$config$speed_breeding) 2L else 4L
     state$config$minimum_mother_age<-4L
     state$config$first_parental_use<-years
     state$config$ocs_backend<-'optiSel_adultcap32_v12'
     state$config$ocs_max_contribution<-job$config$ocs_max_contribution
     state$config$training_stage<-job$config$training_stage
     state$config$training_cohort_stage<-job$config$training_cohort_stage
     state$config$gs_solver<-job$config$gs_solver
     state$config$design_version<-job$config$design_version
     state$prediction_pt_model<-NULL;state$prediction_at_model<-NULL
     options(almond.decision_audit_dir=file.path(job$directory,'training_audit',scenario))
     state$prediction_model<-NULL;state$model_method<-NULL
     # Re-select at the end of the new scenario's year; do not import the
     # pending baseline choice or another scenario's genomic model.
     state$pending<-NULL
     if(method=='Pedigree') {
       if(!requireNamespace('asreml',quietly=TRUE)) stop('ASReml is unavailable.')
       asreml::asreml.options(workspace=job$config$asreml_workspace,trace=FALSE)
     }
     install_timers()
     cat('Scenario',scenario,'\n')
     for(k in seq_len(job$config$future)) {
       state$year<-job$config$burnin+k
       progress('running',scenario,state$year)
       records[[k]]<<-redesign_year(state,method,years,job$rep)
       records[[k]]$scenario<<-scenario
       completed_year<<-state$year
       partial<-redesign_bind(records)
       saveRDS(partial,file.path(job$directory,paste0('results_partial_',scenario,'.rds')))
       if(isTRUE(state$config$ocs_enabled)||state$config$speed_breeding) {
         saveRDS(list(solver=state$speed_ocs_solver,crosses=state$speed_cross_trace),
           file.path(job$directory,sprintf('allocation_%s_year_%03d.rds',scenario,state$year)))
       }
       cat('Completed future year',k,'of',job$config$future,'\n')
     }
     result<-redesign_bind(records)
     stopifnot(nrow(result)==job$config$future,!anyDuplicated(result[c('rep','scenario','year')]))
     saveRDS(data.frame(id=seq_len(SP$lastId),dam=SP$pedigree[,1],sire=SP$pedigree[,2]),
       file.path(job$directory,paste0('pedigree_',scenario,'.rds')))
     saveRDS(result,file.path(job$directory,paste0('results_',scenario,'.rds')))
     result
   })
   statuses[[scenario]]<-data.frame(rep=job$rep,phase=job$phase,scenario=scenario,
     status=if(outcome$ok) 'complete' else 'failed',completed_years=completed_year-job$config$burnin,
     error=if(outcome$ok) NA_character_ else outcome$error,
     elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),
     stringsAsFactors=FALSE)
   saveRDS(redesign_bind(statuses),file.path(job$directory,paste0('scenario_status_',job$phase,'.rds')))
   if(outcome$ok) {
     outputs[[scenario]]<-outcome$result
     progress('scenario_complete',scenario,completed_year)
   } else {
     progress('scenario_failed',scenario,completed_year+1L,outcome$error)
     writeLines(outcome$error,file.path(job$directory,paste0('error_',scenario,'.txt')))
     cat('SCENARIO FAILED; continuing other schemes:',scenario,outcome$error,'\n')
   }
   saveRDS(redesign_bind(outputs),file.path(job$directory,paste0('results_',job$phase,'.rds')))
 }
 manifest<-redesign_bind(statuses)
 failed<-if(is.null(manifest)) FALSE else any(manifest$status=='failed')
 result<-list(ok=TRUE,rep=job$rep,result=redesign_bind(outputs),
   scenario_status=manifest,has_failures=failed)
 progress(if(failed) 'finished_with_failures' else 'complete','all')
 saveRDS(result[names(result)!='result'],file.path(job$directory,paste0('worker_status_',job$phase,'.rds')))
 result
}

