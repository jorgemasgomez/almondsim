# Compatibility entry point for the active revised engine.
run_parallel_replicates <- function(main_dir,config=redesign_config(),scenarios=redesign_all_scenarios(main_dir)) {
 run_redesigned_schemes(main_dir,config,scenarios)
}
