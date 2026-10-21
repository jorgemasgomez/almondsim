# Engine loader used by run_scheme.R, run_experiment.R and scheme-specific 00RUNME.R.
# Scientific logic is defined once in the named R/ modules; no simulations start here.
library(AlphaSimR)
local({
  # source() exposes ofile, whereas sys.source() exposes file in its frame.
  source_files <- unlist(lapply(sys.frames(), function(frame) {
    candidates <- list(frame$ofile, frame$file)
    Filter(function(path) is.character(path) && length(path)==1L &&
      !is.na(path) && basename(path)=='redesigned_schemes.R' && file.exists(path), candidates)
  }), use.names=FALSE)
  if (!length(source_files)) stop('Source this saved file to locate the pipeline modules.')
  pipeline_root <- dirname(normalizePath(tail(source_files,1)[[1]], winslash='/', mustWork=TRUE))
  target <- parent.env(environment())
  modules <- c(
    'simulation_parameters.R',
    'simulate_phenotypes.R',
    'build_training_population.R',
    'fit_prediction_model.R',
    'predict_genetic_values.R',
    'apply_selection.R',
    'select_parents.R',
    'record_results.R',
    'advance_year.R',
    'initialize_pipeline.R',
    'runtime_diagnostics.R',
    'run_replicate.R',
    'run_experiment.R',
    'freeze_pipeline_source.R',
    'run_interface.R'
  )
  for (module in modules) sys.source(file.path(pipeline_root,'R',module), envir=target)
})
