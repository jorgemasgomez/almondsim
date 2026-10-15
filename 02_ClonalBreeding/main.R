# Explicitly running this file starts a fresh experiment; sourcing the engine does not.
.libPaths(c(Sys.getenv('ALMONDSIM_R_LIBRARY','C:/Users/franc/AppData/Local/R/win-library/4.4'),.libPaths()))
files<-Filter(Negate(is.null),lapply(sys.frames(),function(f) f$ofile))
flag<-grep('^--file=',commandArgs(),value=TRUE)
script<-if(length(files)) tail(files,1)[[1]] else if(length(flag)) sub('^--file=','',flag[1]) else NA_character_
if(is.na(script)) stop('Run this saved script so its directory can be resolved.')
root<-dirname(normalizePath(script,winslash='/',mustWork=TRUE))
source(file.path(root,'redesigned_schemes.R'))
config<-redesign_config();scenarios<-redesign_all_scenarios(root)
print(scenarios[c('folder_name','first_parental_use','phenotyping_years','early_selection','ocs_backend')])
run<-run_redesigned_schemes(root,config,scenarios)
source(file.path(root,'plot_results.R'))
plot_experiment(run$run_dir)
