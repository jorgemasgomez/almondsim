# Copy executable sources before starting workers so a running experiment never
# depends on later edits to the working checkout. Output data remain external.
freeze_pipeline_source <- function(original_root, destination) {
  dir.create(destination, recursive=TRUE, showWarnings=FALSE)
  scripts <- list.files(original_root, pattern='[.]R$', full.names=TRUE)
  if (!all(file.copy(scripts, destination))) stop('Could not freeze top-level R scripts.')
  module_source <- file.path(original_root, 'R')
  module_destination <- file.path(destination, 'R')
  dir.create(module_destination, recursive=TRUE, showWarnings=FALSE)
  modules <- list.files(module_source, pattern='[.]R$', full.names=TRUE)
  if (!length(modules) || !all(file.copy(modules, module_destination))) {
    stop('Could not freeze the breeding pipeline modules.')
  }
  # Relative paths make the manifest usable after moving a completed run.
  burnin_destination <- file.path(destination,'00_Burn_in')
  dir.create(burnin_destination,recursive=TRUE,showWarnings=FALSE)
  if (!file.copy(file.path(original_root,'00_Burn_in','GlobalParameters.R'),
                 burnin_destination,overwrite=TRUE)) stop('Could not freeze burn-in parameters.')
  files <- c(basename(scripts), file.path('R', basename(modules)),
             file.path('00_Burn_in','GlobalParameters.R'))
  checksums <- unname(tools::md5sum(file.path(destination, files)))
  write.csv(data.frame(file=files, md5=checksums),
            file.path(destination, 'source_manifest.csv'), row.names=FALSE)
  invisible(destination)
}
