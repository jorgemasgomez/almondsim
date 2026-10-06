# Temporary two-replicate comparison using the historical haplotypes.
# Supply founder count explicitly: Rscript run_ocs_legacy_pilot.R 20
.libPaths(c('C:/Users/franc/AppData/Local/R/win-library/4.4',.libPaths()))
root <- 'C:/Users/franc/Desktop/Proyectos/Almond_sim/almondsim/02_ClonalBreeding'
source(file.path(root,'redesigned_schemes.R'))
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=1L) stop('Specify the number of diploid founders explicitly.')
founder_count <- as.integer(arguments[1])
stopifnot(founder_count %in% c(10L,20L,30L))
config <- redesign_config()
config$founder_source <- 'legacy';config$legacy_founders <- founder_count
config$reps <- 2L;config$workers <- 2L;config$threads <- 2L
scenarios <- redesign_scenarios(root)
names <- c('GS_Y2_Early1500',paste0('GS_Y2_Early1500_OCS_',c('Diversity','Balanced','Gain')))
scenarios <- scenarios[scenarios$scenario_id %in% names,,drop=FALSE]
stopifnot(nrow(scenarios)==4L)
# Run these two pilot replicates in the main R process. This avoids a
# local socket-worker startup issue and shares fewer resources with stdpopsim.
config$workers <- 1L
run <- file.path(redesign_data_dir(root,config),'redesigned_runs',paste0('ocs_legacy_pilot_',format(Sys.time(),'%Y%m%d_%H%M%S')))
dir.create(run,recursive=TRUE)
writeLines(run,file.path(root,'ocs_legacy_pilot_latest.txt'))
saveRDS(list(config=config,scenarios=scenarios),file.path(run,'configuration.rds'))
RNGkind("L'Ecuyer-CMRG");set.seed(config$seed)
streams <- list(.Random.seed,parallel::nextRNGStream(.Random.seed))
results <- vector('list',2L)
for(rep in 1:2) {
  jsonlite::write_json(list(state='running',rep=rep,replicates=2,run=run),
                       file.path(run,'pilot_status.json'),pretty=TRUE,auto_unbox=TRUE)
  job <- list(rep=rep,root=root,config=config,scenarios=scenarios,
              methods='GS',phase='parallel',seed=streams[[rep]],libraries=.libPaths(),
              directory=file.path(run,sprintf('rep_%03d',rep)))
  cat('Starting pilot replicate',rep,'in',run,'\n');flush.console()
  results[[rep]] <- redesign_worker(job)
  if(!results[[rep]]$ok) {
    jsonlite::write_json(list(state='failed',rep=rep,error=results[[rep]]$error),
                         file.path(run,'pilot_status.json'),pretty=TRUE,auto_unbox=TRUE)
    stop(results[[rep]]$error)
  }
}
combined <- do.call(rbind,lapply(results,`[[`,'result'))
saveRDS(combined,file.path(run,'results_all_schemes.rds'))
write.table(combined,file.path(run,'results_all_schemes.txt'),sep='\t',row.names=FALSE,quote=FALSE)
saveRDS(list(complete=TRUE,replicates=2L,scenarios=scenarios),file.path(run,'status.rds'))
jsonlite::write_json(list(state='complete',replicates=2,run=run),
                     file.path(run,'pilot_status.json'),pretty=TRUE,auto_unbox=TRUE)
cat('PILOT_COMPLETE',run,'\n')
