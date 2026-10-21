# Record annual outcomes by population stage and combine records with compatible columns.
# Public function names are retained for existing launchers and saved workflows.

# Record population sizes, genetic means, genetic variances and selection accuracies for the current year.
redesign_record <- function(state,scenario,rep,accuracy,entrants,parent_stage) {
 row <- list(year=state$year,rep=rep,scenario=scenario,accSel=accuracy,
             accParents=state$parent_accuracy,newParents=entrants,parentSelectionStage=parent_stage)
 for(stage in c('Parents',state$config$stages)) {
  pop <- if(stage=='Parents') state$parents else state$pop[[stage]]
  row[[paste0('n_',stage)]] <- if(is.null(pop)) 0L else nInd(pop)
  row[[paste0('meanG_',stage)]] <- if(is.null(pop)) NA_real_ else as.numeric(meanG(pop))[1L]
  row[[paste0('varG_',stage)]] <- if(is.null(pop)) NA_real_ else as.numeric(varG(pop))[1L]
 }
 row$meanG <- row$meanG_Seedlings;row$varG <- row$varG_Seedlings
 haplo <- pullMarkerHaplo(state$pop$Seedlings,markers=locuscompt_position_vector,haplo='all',simParam=SP)
 row$allelesSI <- length(unique(apply(haplo,1,paste,collapse='_')))
 as.data.frame(row,stringsAsFactors=FALSE)
}

# Combine annual/status frames while preserving all columns when some records lack a field.
redesign_bind <- function(frames) {
 frames<-Filter(function(x) !is.null(x)&&nrow(x)>0L,frames)
 if(!length(frames)) return(NULL)
 columns<-unique(unlist(lapply(frames,names)))
 do.call(rbind,lapply(frames,function(x) {
   for(column in setdiff(columns,names(x))) x[[column]]<-NA
   x[columns]
 }))
}

