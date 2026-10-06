# Run only this scenario; use project main.R for all 25 scenarios.
.libPaths(c('C:/Users/franc/AppData/Local/R/win-library/4.4',.libPaths()))
script_files <- Filter(Negate(is.null),lapply(sys.frames(),function(frame) frame$ofile))
script_flag <- grep('^--file=',commandArgs(),value=TRUE)
script_file <- if(length(script_files)) tail(script_files,1)[[1]] else if(length(script_flag)) sub('^--file=','',script_flag[1]) else NA_character_
if(is.na(script_file)) stop('Run or source this saved script so its project directory can be resolved.')
main_dir <- dirname(dirname(normalizePath(script_file,winslash='/',mustWork=TRUE)))
source(file.path(main_dir,'redesigned_schemes.R'))
config <- redesign_config()
scenarios <- redesign_scenarios(main_dir)
scenarios <- scenarios[scenarios$folder_name=='16_GS_2_Early1500',,drop=FALSE]
stopifnot(nrow(scenarios)==1L)
parallel_run <- run_redesigned_schemes(main_dir,config,scenarios)
