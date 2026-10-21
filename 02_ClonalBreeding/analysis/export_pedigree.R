# Explicit exporter; no file is read or generated merely by sourcing this file.
export_pedigree <- function(run,rep,scenario) {
 file<-file.path(run,sprintf('rep_%03d',rep),paste0('pedigree_',scenario,'.rds'))
 if(!file.exists(file)) stop('No saved completed pedigree for this scenario.')
 pedigree<-readRDS(file)
 destination<-file.path(run,sprintf('rep_%03d',rep),paste0('pedigree_',scenario,'.csv'))
 write.csv(pedigree,destination,row.names=FALSE)
 invisible(destination)
}
