# Archive evaluated cohorts and construct decision-specific training sets. Do not replace an original PT cohort by its selected AT survivors.
# Function names describe each operation; all launchers use this shared implementation.

# Merge population versions by ID, preferring the most recently supplied version.
merge_latest_populations <- function(pops) {
  pops <- Filter(Negate(is.null),pops)
  combined <- mergePops(pops)
  combined[which(!duplicated(as.character(combined@id),fromLast=TRUE))]
}

# Save complete evaluated cohorts before selection, retaining cohort identity and record year.
archive_evaluated_cohorts <- function(state) {
  if (is.null(state$stage_history)) {
    stop('This burn-in predates stage-specific training archives. Generate a new burn-in; do not reuse historical caches.')
  }
  for (stage in c('PT2','AT3','ET3')) {
    pop <- state$pop[[stage]]
    if (is.null(pop)) next
    entries <- state$stage_history[[stage]]
    # Replace, rather than duplicate, an archive if the same year is recorded twice.
    entries[[as.character(state$year)]] <- pop
    lag <- match('AT3',state$config$stages)-match('PT2',state$config$stages)
    keep <- state$config$training_cohorts+if(stage=='PT2') lag else 0L
    state$stage_history[[stage]] <- tail(entries,keep)
  }
  invisible(state)
}

# Build the PT-only training population; retained as an explicit implementation, not a source-order alias.
build_pt_training_population <- function(state,candidates,anchor=NULL) {
  stage <- if(is.null(anchor)) state$training_cohort_stage else anchor
  if(is.null(stage)) stop('Missing cohort training anchor.')
  cohorts <- tail(state$stage_history[[stage]],state$config$training_cohorts)
  if(length(cohorts)!=state$config$training_cohorts) stop('Training requires three complete archived cohorts.')
  lag <- match('AT3',state$config$stages)-match('PT2',state$config$stages)
  full <- lapply(names(cohorts),function(year) {
    if(stage=='AT3') {
      baseline <- state$stage_history$PT2[[as.character(as.integer(year)-lag)]]
      survivors <- cohorts[[year]]
      if(is.null(baseline)) stop('Missing matching full PT2 archive. Generate a new burn-in with six PT2 archives; no transition or substitute cohort is allowed.')
      if(nInd(baseline)!=state$config$seedlings_keep || nInd(survivors)!=state$config$act ||
         anyDuplicated(baseline@id) || anyDuplicated(survivors@id) ||
         !all(survivors@id %in% baseline@id)) stop('PT2/AT3 cohort identity or size mismatch.')
      # AT survivors validate cohort identity only. ALL 500 original PT records remain.
      pop <- baseline
      labels <- rep('PT',nInd(pop))
    } else {
      if(stage!='PT2') stop('Unsupported training stage.')
      pop <- cohorts[[year]]; labels <- rep('PT',nInd(pop))
      if(nInd(pop)!=state$config$seedlings_keep) stop('PT2 training cohort is not complete.')
    }
    list(pop=pop,stage=labels,cohort=rep(year,nInd(pop)))
  })
  train <- merge_latest_populations(lapply(full,`[[`,'pop'))
  ids <- unlist(lapply(full,function(x)as.character(x$pop@id)),use.names=FALSE)
  if(anyDuplicated(ids) || nInd(train)!=state$config$training_cohorts*state$config$seedlings_keep) stop('Training requires 1500 unique individuals.')
  positions <- match(as.character(train@id),ids)
  train@misc$training_record_stage <- unlist(lapply(full,`[[`,'stage'),use.names=FALSE)[positions]
  train@misc$training_record_cohort <- unlist(lapply(full,`[[`,'cohort'),use.names=FALSE)[positions]
  counts <- table(train@misc$training_record_stage)
  if(!identical(as.integer(counts['PT']),1500L) || any(train@misc$training_record_stage!='PT')) stop('Training must contain only 1500 original PT records.')
  train
}

# Identify schemes whose AT/parent decisions require calibrated repeated PT and AT records.
uses_pt_at_records <- function(state) isTRUE(state$config$first_parental_use==9L)

# Dispatch training construction according to the decision and repeated-record design.
build_training_population <- function(state,candidates,anchor=NULL) {
  stage <- if(is.null(anchor)) state$training_cohort_stage else anchor
  train <- build_pt_training_population(state,candidates,anchor=stage)
  extra <- data.frame(id=character(),pheno=numeric(),record_year=integer(),stage=character(),cohort=character())
  if(uses_pt_at_records(state) && stage=='AT3') {
    cohorts <- tail(state$stage_history$AT3,state$config$training_cohorts)
    extra <- do.call(rbind,lapply(names(cohorts),function(year) {
      p <- cohorts[[year]]
      data.frame(id=as.character(p@id),pheno=as.numeric(p@pheno),record_year=as.integer(p@fixEff),stage='AT',cohort=year)
    }))
    stopifnot(nrow(extra)==75L,!anyDuplicated(extra$id),all(extra$id %in% train@id))
    idx <- match(extra$id,as.character(train@id))
    stopifnot(all(extra$cohort==train@misc$training_record_cohort[idx]),all(extra$record_year-train@fixEff[idx]==3L))
  }
  train@misc$extra_AT_records <- extra
  train
}

# Construct the repeated-record table while checking the required cohort and observation counts.
build_pt_at_records <- function(train) {
  pt <- data.frame(id=as.character(train@id),pheno=as.numeric(train@pheno),record_year=as.integer(train@fixEff),stage='PT',cohort=as.character(train@misc$training_record_cohort))
  records <- rbind(pt,train@misc$extra_AT_records)
  stopifnot(sum(records$stage=='PT')==1500L,length(unique(records$id))==1500L,!anyDuplicated(records[c('id','stage')]))
  if(nrow(train@misc$extra_AT_records))stopifnot(nrow(records)==1575L,sum(table(records$id)==2L)==75L,sum(records$stage=='AT')==75L)
  records
}

# Expand prediction candidates onto PT/AT records, retaining missing responses for unobserved candidates.
build_pedigree_record_frame <- function(records,prediction) {
  ids <- sort(unique(as.character(prediction@id)))
  frame <- expand.grid(Stage=c('PT','AT'),Ind=ids,stringsAsFactors=FALSE)
  frame$Ind <- factor(frame$Ind,levels=ids);frame$Stage <- factor(frame$Stage,levels=c('PT','AT'))
  frame <- frame[order(frame$Ind,frame$Stage),,drop=FALSE]
  key <- function(id,stage)paste(id,stage,sep=':')
  positions <- match(key(as.character(frame$Ind),as.character(frame$Stage)),key(records$id,records$stage))
  frame$Pheno <- records$pheno[positions]
  years <- records$record_year[positions]
  years[is.na(years)] <- min(records$record_year)
  frame$Year <- factor(years)
  stopifnot(sum(!is.na(frame$Pheno))==1575L,all(table(frame$Ind)==2L))
  frame
}

