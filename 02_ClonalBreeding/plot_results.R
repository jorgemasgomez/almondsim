# Default manuscript plotting entry point; sourcing never launches a simulation.
plot_experiment <- local({
 sources<-Filter(Negate(is.null),lapply(sys.frames(),function(f)f$ofile))
 module<-if(length(sources))file.path(dirname(tail(sources,1)[[1]]),'plot_manuscript.R')else NA_character_
 function(run,with_intervals=TRUE,...) {
  file<-getOption('almond.manuscript_plot_source')
  if(is.null(file)) {
   candidates<-c(module,file.path(run,'code','plot_manuscript.R'),
    file.path(getOption('almond.input_dir',''),'plot_manuscript.R'),
    file.path('02_ClonalBreeding','plot_manuscript.R'),'plot_manuscript.R')
   existing<-candidates[!is.na(candidates)&file.exists(candidates)]
   if(!length(existing))stop('Cannot locate plot_manuscript.R; source it explicitly and call plot_manuscript().')
   file<-existing[1L]
  }
  env<-new.env(parent=environment());sys.source(file,envir=env)
  env$plot_manuscript(run,with_intervals=with_intervals,...)
 }
})
