# Freeze source files and coordinate replicate workers. Pedigree processing remains serial within the existing execution policy.
# Public function names are retained for existing launchers and saved workflows.

# Start a new experiment only when explicitly called. Sourcing the module does not start a run.
run_redesigned_schemes <- function(root,config=redesign_config(),scenarios=redesign_scenarios()) {
 experiment_start<-proc.time()[['elapsed']]
 original_root<-normalizePath(root,winslash='/',mustWork=TRUE)
 if(!identical(config$founder_source,'bank150')) stop('The revised design requires bank150.')
 redesign_founder_bank(config$founder_bank)
 config$data_dir<-redesign_data_dir(original_root,config)
 if(is.null(config$asreml_workspace)) config$asreml_workspace<-'16gb'
 parent<-file.path(config$data_dir,'manuscript_runs')
 dir.create(parent,recursive=TRUE,showWarnings=FALSE)
 run<-tempfile(paste0('PT_AT_ET_',config$reps,'reps_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),tmpdir=parent)
 dir.create(run,recursive=TRUE)
 # Freeze source so later manuscript/code edits cannot change a running run.
 root<-file.path(run,'code');dir.create(root)
 dir.create(file.path(root,'00_Burn_in'))
 # Freeze every executable module as well as the compatibility entry points.
 freeze_pipeline_source(original_root,root)
 writeLines(run,file.path(config$data_dir,'latest_run.txt'))
 RNGkind("L'Ecuyer-CMRG");set.seed(config$seed);stream<-.Random.seed
 jobs<-vector('list',config$reps)
 for(rep in seq_len(config$reps)) {
   jobs[[rep]]<-list(rep=rep,root=root,config=config,scenarios=scenarios,
     methods=c('Pheno','GS'),phase='parallel',seed=stream,libraries=.libPaths(),
     directory=file.path(run,sprintf('rep_%03d',rep)))
   dir.create(jobs[[rep]]$directory)
   saveRDS(jobs[[rep]],file.path(jobs[[rep]]$directory,'job_parallel.rds'))
   stream<-parallel::nextRNGStream(stream)
 }
 saveRDS(list(config=config,scenarios=scenarios,source_root=original_root),file.path(run,'configuration.rds'))
 write.csv(scenarios,file.path(run,'scenarios.csv'),row.names=FALSE)
 jsonlite::write_json(list(status='running',phase='parallel',replicates=config$reps,
   schemes=nrow(scenarios),started=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')),
   file.path(run,'experiment_progress.json'),auto_unbox=TRUE,pretty=TRUE)
 cluster<-parallel::makeCluster(min(config$workers,config$reps),outfile=file.path(run,'cluster.log'))
 on.exit(parallel::stopCluster(cluster),add=TRUE)
 results<-parallel::parLapplyLB(cluster,jobs,function(job) {
   .libPaths(job$libraries);source(file.path(job$root,'redesigned_schemes.R'));redesign_worker(job)
 })
 # Scheme failures do not stop unrelated schemes, replicates or pedigree.
 jsonlite::write_json(list(status='running',phase='pedigree',replicates=config$reps,
   schemes=nrow(scenarios)),file.path(run,'experiment_progress.json'),auto_unbox=TRUE,pretty=TRUE)
 pedigree<-vector('list',config$reps)
 for(rep in seq_len(config$reps)) {
   job<-jobs[[rep]];job$methods<-'Pedigree';job$phase<-'pedigree'
   pedigree[[rep]]<-parallel::parLapply(cluster[1],list(job),function(job) {
     .libPaths(job$libraries);source(file.path(job$root,'redesigned_schemes.R'));redesign_worker(job)
   })[[1]]
 }
 outcomes<-c(results,pedigree)
 combined<-redesign_bind(lapply(outcomes,function(x) x$result))
 statuses<-redesign_bind(lapply(outcomes,function(x) x$scenario_status))
 fatal<-Filter(function(x) !isTRUE(x$ok),outcomes)
 complete<-!length(fatal)&&(nrow(scenarios)==0L || (!is.null(statuses)&&
   all(statuses$status=='complete')&&nrow(statuses)==config$reps*nrow(scenarios)))
 expected_rows<-config$reps*nrow(scenarios)*config$future
 if(complete && expected_rows>0L) stopifnot(nrow(combined)==expected_rows,!anyDuplicated(combined[c('rep','scenario','year')]))
 saveRDS(combined,file.path(run,'results_all_schemes.rds'))
 saveRDS(statuses,file.path(run,'scenario_status_all.rds'))
 saveRDS(list(complete=complete,finished=TRUE,replicates=config$reps,schemes=nrow(scenarios),
   expected_rows=expected_rows,verified_rows=if(is.null(combined)) 0L else nrow(combined),
   failed_scenarios=if(is.null(statuses)) NULL else statuses[statuses$status=='failed',],
   initialization_failures=lapply(fatal,function(x) x[c('rep','error')])),file.path(run,'status.rds'))
 saveRDS(list(elapsed_seconds=proc.time()[['elapsed']]-experiment_start,completed=Sys.time()),
   file.path(run,'total_runtime.rds'))
 jsonlite::write_json(list(status=if(complete) 'complete' else 'finished_with_failures',
   replicates=config$reps,schemes=nrow(scenarios),expected_rows=expected_rows,
   verified_rows=if(is.null(combined)) 0L else nrow(combined)),
   file.path(run,'experiment_progress.json'),auto_unbox=TRUE,pretty=TRUE)
 message('Experiment finished: ',run,'; all scenarios complete: ',complete)
 invisible(list(run_dir=run,data=combined,statuses=statuses))
}

