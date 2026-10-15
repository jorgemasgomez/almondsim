# Explicitly running this saved file launches only GS_SpeedBreeding_Y2_Early1500_OCS_Balanced.
.libPaths(c(Sys.getenv('ALMONDSIM_R_LIBRARY','C:/Users/franc/AppData/Local/R/win-library/4.4'),.libPaths()))
files<-Filter(Negate(is.null),lapply(sys.frames(),function(f) f$ofile))
flag<-grep('^--file=',commandArgs(),value=TRUE)
script<-if(length(files)) tail(files,1)[[1]] else if(length(flag)) sub('^--file=','',flag[1]) else NA_character_
if(is.na(script)) stop('Run the saved launcher to resolve its folder.')
root<-dirname(dirname(normalizePath(script,winslash='/',mustWork=TRUE)))
source(file.path(root,'redesigned_schemes.R'))
config<-redesign_config();scenarios<-redesign_all_scenarios(root)
scenarios<-scenarios[scenarios$scenario_id=='GS_SpeedBreeding_Y2_Early1500_OCS_Balanced',,drop=FALSE]
stopifnot(nrow(scenarios)==1L)
run<-run_redesigned_schemes(root,config,scenarios)
source(file.path(root,'plot_results.R'));plot_experiment(run$run_dir)
