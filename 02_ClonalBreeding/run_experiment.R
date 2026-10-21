# Human-facing launcher for a selected set of schemes. No arguments prints help.
# This shares the engine and annual cycle with run_scheme.R and legacy launchers.
.libPaths(c(Sys.getenv('ALMONDSIM_R_LIBRARY',
                      'C:/Users/franc/AppData/Local/R/win-library/4.4'), .libPaths()))
.pipeline_files <- unlist(lapply(sys.frames(), function(frame) {
  Filter(function(path) is.character(path)&&length(path)==1L&&!is.na(path)&&
    basename(path)=='run_experiment.R'&&file.exists(path),list(frame$ofile,frame$file))
}),use.names=FALSE)
.pipeline_flag <- grep('^--file=', commandArgs(), value=TRUE)
.pipeline_script <- if(length(.pipeline_files)) tail(.pipeline_files,1)[[1]] else
  if(length(.pipeline_flag)) sub('^--file=','',.pipeline_flag[1]) else NA_character_
if(is.na(.pipeline_script)) stop('Run or source this saved script to resolve its directory.')
.pipeline_root <- dirname(normalizePath(.pipeline_script,winslash='/',mustWork=TRUE))
source(file.path(.pipeline_root,'redesigned_schemes.R'),local=environment())
if (sys.nframe()==0L) pipeline_cli(.pipeline_root, single=FALSE)
