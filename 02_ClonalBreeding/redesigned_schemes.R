# Revised manuscript design: first parental use in years 9, 6, 4 or 2.
# PT (4-5), AT (6-8), ET (9-11); GS/Ped9 repeated calibrated PT+AT.
library(AlphaSimR)

redesign_config <- function() {
 list(design_version='PT_AT_ET_optiSel_adultcap32_v12',reps=10L,workers=4L,threads=6L,
  seed=20261005L,asreml_workspace='16gb',burnin=40L,future=40L,fill_cycles=14L,
  parents=32L,crosses=20L,progeny=25L,act=25L,ect=8L,entrants=8L,
  training_cohorts=3L,gs_solver='AlphaSimR::RRBLUP',early_selection=FALSE,seedlings_keep=500L,
  ocs_alpha=NA_real_,ocs_target_Ne=NA_real_,ocs_max_contribution=.10,ocs_family_sizes='variable',ocs_max_parents=NA_integer_,ocs_backend='optiSel_adultcap32_v12',
  speed_breeding=FALSE,minimum_mother_age=4L,minimum_pollen_age=4L,
  data_dir=Sys.getenv('ALMONDSIM_DATA_DIR',unset=''),founder_source='bank150',
  founder_bank=Sys.getenv('ALMONDSIM_FOUNDER_BANK',unset=
    'C:/Users/franc/Desktop/Proyectos/Almond_sim/Haplotipos_2026-10-05_60000'),
  stages=c('F1','Seedlings','Juvenile2','Juvenile3','PT1','PT2',
           'AT1','AT2','AT3','ET1','ET2','ET3'),
  h2=c(PT1=.2,PT2=.2,AT1=.4,AT2=.4,AT3=.4,ET1=.5,ET2=.5,ET3=.5,Parents=.5),
  parent_stage=c(`9`='AT3',`6`='PT2',`4`='PT1',`2`='Juvenile2'),
  training_stage=c(`9`='PT2',`6`='PT2',`4`='PT2',`2`='PT2'),
  training_cohort_stage=c(`9`='AT3',`6`='PT2',`4`='PT2',`2`='PT2'))
}

redesign_scenarios <- function(root=NULL) {
 make <- function(method,first,early=FALSE,renewal=25,alpha=NA_real_,profile=NA_character_) {
  speed <- first==2L
  suffix <- if(!is.na(alpha)) paste0('_OCS_',profile) else if(renewal!=25) paste0('_R',renewal) else ''
  stem <- paste0(method,if(speed) '_SpeedBreeding' else '', '_Y',first,
                 if(early) '_Early1500' else '',suffix)
  data.frame(method=method,first_parental_use=as.integer(first),
   phenotyping_years=if(first==9L) 5L else if(first==6L) 2L else 0L,
   parent_selection_stage=redesign_config()$parent_stage[[as.character(first)]],
   training_stage=redesign_config()$training_stage[[as.character(first)]],
   training_cohort_stage=redesign_config()$training_cohort_stage[[as.character(first)]],
   selection_timing=if(first>=6L) 'end_of_previous_year' else 'before_own_phenotyping',
   early_selection=early,progeny_per_cross=if(early) 75L else 25L,
   seedlings_keep=500L,new_parents_per_year=if(is.na(alpha)) as.integer(32*renewal/100) else NA_integer_,
   renewal_percent=if(is.na(alpha)) renewal else NA_real_,ocs_alpha=alpha,
   ocs_profile=profile,ocs_target_Ne=if(is.na(alpha)) NA_real_ else unname(c(Diversity=30,Balanced=10,Gain=5)[profile]),ocs_backend=if(is.na(alpha)) NA_character_ else 'optiSel_adultcap32_v12',
   family_size_policy=if(is.na(alpha)) 'equal' else 'variable',
   target_offspring=if(early) 1500L else 500L,
   ocs_max_contribution=if(is.na(alpha)) NA_real_ else .10,
   ocs_max_parents=NA_integer_,
   speed_breeding=speed,scenario_id=stem,stringsAsFactors=FALSE)
 }
 rows <- list()
 for(first in c(9L,6L)) for(method in c('Pheno','Pedigree','GS')) rows[[length(rows)+1L]] <- make(method,first)
 for(method in c('Pedigree','GS')) rows[[length(rows)+1L]] <- make(method,4L)
 rows[[length(rows)+1L]] <- make('GS',2L)
 for(method in c('Pedigree','GS')) rows[[length(rows)+1L]] <- make(method,6L,TRUE)
 rows[[length(rows)+1L]] <- make('GS',2L,TRUE)
 for(rate in c(50,75,100)) rows[[length(rows)+1L]] <- make('GS',6L,FALSE,rate)
 for(first in c(6L,2L)) for(i in 1:3) rows[[length(rows)+1L]] <-
   make('GS',first,TRUE,alpha=c(.10,.20,.40)[i],profile=c('Diversity','Balanced','Gain')[i])
 for(i in 1:3) rows[[length(rows)+1L]] <- make('GS',4L,TRUE,alpha=c(.10,.20,.40)[i],profile=c('Diversity','Balanced','Gain')[i])
 result <- do.call(rbind,rows)
 result$folder_name <- sprintf('%02d_%s',seq_len(nrow(result)),result$scenario_id)
 if(!is.null(root)) {
  for(i in seq_len(nrow(result))) {
   settings <- new.env(parent=baseenv())
   sys.source(file.path(root,result$folder_name[i],'GlobalParameters.R'),envir=settings)
   for(field in c('selection_method','first_parental_use','parent_selection_stage',
                 'training_stage','early_selection','progeny_per_cross','ocs_alpha','speed_breeding')) {
    key <- if(field=='selection_method') 'method' else field
    if(!identical(settings[[field]],result[[key]][i])) stop('Scenario settings differ from the catalogue: ',result$folder_name[i],' / ',field)
   }
  }
  result$folder <- file.path(root,result$folder_name)
 }
 result
}

# Nested means: covariance between successive aggregate errors equals the
# smaller error variance. Equal target h2 retains the previous aggregate.
# Initial h2 is calibrated against founder additive variance. Error variance
# remains fixed thereafter; realized h2 can decline with genetic variance.
# Correlated aggregate residuals with explicitly calibrated marginal variances.
# Annual GxE is generated separately; raw annual residual records are not stored.
redesign_error_step <- function(previous,old_variance,new_variance) {
  if(old_variance<=1e-12) return(rnorm(length(previous),sd=sqrt(new_variance)))
  covariance <- min(old_variance,new_variance)
  ratio <- covariance/old_variance
  ratio*previous+rnorm(length(previous),sd=sqrt(max(0,new_variance-covariance^2/old_variance)))
}

redesign_environment_schedule <- function(config,seed) {
  had_seed <- exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) old_seed <- get('.Random.seed',envir=.GlobalEnv)
  on.exit(if(had_seed) assign('.Random.seed',old_seed,envir=.GlobalEnv) else
    if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv))
  set.seed(seed)
  years <- (-config$fill_cycles):(config$burnin+config$future)
  setNames(runif(length(years)),as.character(years))
}

redesign_measure <- function(pop,stage,state) {
  if(is.null(pop)) return(NULL)
  if(is.null(state$environment_p) || is.null(state$annual_evaluations))
    stop('Missing annual environment state. Generate a new burn-in.')
  probability <- state$environment_p[[as.character(state$year)]]
  if(is.null(probability) || !is.finite(probability)) stop('Annual environment is missing for this year.')
  trait <- SP$traits[[1]]
  if(!.hasSlot(trait,'envVar') || !length(pop@gxe)) stop('Annual evaluation requires the original GxE trait.')
  annual_z <- qnorm(probability,sd=sqrt(trait@envVar))
  residuals <- environment <- numeric(nInd(pop))
  ids <- as.character(pop@id)
  for(i in seq_along(ids)) {
    history <- state$annual_evaluations[[ids[i]]]
    if(is.null(history)) history <- list(years=integer(),sum_z=0,n=0L)
    if(!state$year %in% history$years) {
      history$years <- c(history$years,state$year)
      history$sum_z <- history$sum_z+annual_z;history$n <- history$n+1L
      state$annual_evaluations[[ids[i]]] <- history
    }
    environment[i] <- pop@gxe[[1]][i]*history$sum_z/history$n
    # Budget includes averaged GxE: do not add it on top of the full target error.
    total_error <- state$reference_variance*(1/state$config$h2[[stage]]-1)
    variance <- total_error-state$reference_gxe_variance/history$n
    if(variance < -1e-8) stop('Target aggregate h2 is incompatible with annual GxE variance.')
    variance <- max(0,variance)
    prior <- state$errors[[ids[i]]]
    if(is.null(prior)) error <- rnorm(1,sd=sqrt(variance)) else
      error <- redesign_error_step(prior$error,prior$variance,variance)
    state$errors[[ids[i]]] <- list(error=error,variance=variance)
    residuals[i] <- error
  }
  pop@pheno <- matrix(as.numeric(gv(pop))+environment+residuals,ncol=1)
  pop@fixEff <- rep(as.integer(state$year),nInd(pop))
  pop
}

redesign_merge_latest <- function(pops) {
  pops <- Filter(Negate(is.null),pops)
  combined <- mergePops(pops)
  combined[which(!duplicated(as.character(combined@id),fromLast=TRUE))]
}

# Archive complete cohorts before the scenario-specific 500 -> 25 cut.
# Separate stage archives allow paired scenarios to share a common burn-in.
redesign_archive_cohorts <- function(state) {
  if (is.null(state$stage_history)) {
    stop('This burn-in predates stage-specific training archives. Generate a new burn-in; do not reuse historical caches.')
  }
  for (stage in setdiff(names(state$config$h2),'Parents')) {
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

redesign_training <- function(state,candidates,anchor=NULL) {
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
  train <- redesign_merge_latest(lapply(full,`[[`,'pop'))
  ids <- unlist(lapply(full,function(x)as.character(x$pop@id)),use.names=FALSE)
  if(anyDuplicated(ids) || nInd(train)!=state$config$training_cohorts*state$config$seedlings_keep) stop('Training requires 1500 unique individuals.')
  positions <- match(as.character(train@id),ids)
  train@misc$training_record_stage <- unlist(lapply(full,`[[`,'stage'),use.names=FALSE)[positions]
  train@misc$training_record_cohort <- unlist(lapply(full,`[[`,'cohort'),use.names=FALSE)[positions]
  counts <- table(train@misc$training_record_stage)
  if(!identical(as.integer(counts['PT']),1500L) || any(train@misc$training_record_stage!='PT')) stop('Training must contain only 1500 original PT records.')
  train
}

# PT-only schemes retain AlphaSimR::RRBLUP. Y9 calibrated models are below.
redesign_native_rrblup <- function(train,state) {
  if(any(train@misc$training_record_stage!='PT')) stop('GS training must contain only original PT records.')
  model <- AlphaSimR::RRBLUP(train,simParam=SP)
  attr(model,'training_diagnostics') <- list(solver='AlphaSimR::RRBLUP',
    version=as.character(packageVersion('AlphaSimR')),record_counts=table(train@misc$training_record_stage),
    cohort_keys=unique(train@misc$training_record_cohort),fixed_year=TRUE)
  model
}

# A decision must use the candidates' own ORIGINAL PT records.
# PT2: current PT cohort and previous two; AT3: original PT records of
# current AT cohort and previous two AT cohorts. Y9 adds repeated AT records below.
redesign_assert_decision_training <- function(pop,model,anchor) {
  if(is.null(pop)) return(invisible(TRUE))
  meta <- attr(model,'decision_training')
  if(is.null(meta) || !identical(meta$anchor,anchor) ||
     !all(as.character(pop@id) %in% meta$ids))
    stop('Decision model excludes own PT records or uses the wrong cohort anchor: ',anchor)
  invisible(TRUE)
}

redesign_fit_for_stage <- function(state,method,anchor) {
  if(method=='Pheno') return(NULL)
  if(!anchor %in% c('PT2','AT3')) stop('Unsupported decision training anchor.')
  old_anchor <- state$training_cohort_stage;old_stage <- state$training_stage
  on.exit({state$training_cohort_stage<-old_anchor;state$training_stage<-old_stage},add=TRUE)
  state$training_cohort_stage <- anchor;state$training_stage <- 'PT2'
  model <- redesign_fit(state,method,state$pop[[anchor]])
  train <- redesign_training(state,NULL,anchor=anchor)
  meta <- list(anchor=anchor,ids=as.character(train@id),
    cohort_keys=unique(train@misc$training_record_cohort),
    record_years=train@fixEff,PT_records=nInd(train),AT_records=0L,
    diagnostics=attr(model,'training_diagnostics'))
  attr(model,'decision_training') <- meta
  audit <- getOption('almond.decision_audit_dir')
  if(!is.null(audit)) {
    dir.create(audit,recursive=TRUE,showWarnings=FALSE)
    last_cohort <- tail(meta$cohort_keys,1)
    saveRDS(meta,file.path(audit,sprintf('fit_year_%03d_%s_last_%s.rds',state$year,anchor,last_cohort)))
  }
  model
}

redesign_trial_cuts <- function(state,method,models) {
 for(entry in list(c('PT2',state$config$act),c('AT3',state$config$ect))) {
  stage <- entry[1];n <- as.integer(entry[2]);pop <- state$pop[[stage]]
  if(!is.null(pop)) {
   model <- if(method=='Pheno') NULL else models[[stage]]
   if(method!='Pheno') redesign_assert_decision_training(pop,model,stage)
   state$pop[[stage]] <- redesign_choose(pop,min(n,nInd(pop)),method,model)
  }
 }
 invisible(state)
}

redesign_fit <- function(state,method,candidates) {
  if (method=='Pheno') return(NULL)
  train <- redesign_training(state,candidates)
  if (method=='GS') return(redesign_native_rrblup(train,state))
  if (identical(getOption('almond.pedigree_solver','asreml'),'validation')) {
    return(redesign_validation_blup(state,train))
  }
  if (!requireNamespace('asreml',quietly=TRUE)) stop('Pedigree requires ASReml.')
  prediction <- redesign_merge_latest(c(list(train,state$parents),state$pop))
  ped <- data.frame(id=as.character(seq_len(SP$lastId)),
                    dam=as.character(SP$pedigree[,1]),sire=as.character(SP$pedigree[,2]))
  ped$dam[ped$dam=='0'] <- NA
  ped$sire[ped$sire=='0'] <- NA
  # Keep all ancestors needed for prediction, including unphenotyped seedlings.
  observed <- match(as.character(prediction@id),as.character(train@id))
  frame <- data.frame(Ind=factor(as.character(prediction@id)),
                      Pheno=as.numeric(train@pheno)[observed],
                      Year=factor(train@fixEff[observed]))
  # Missing response rows still need a valid fixed-effect design for prediction.
  frame$Year[is.na(frame$Year)] <- levels(frame$Year)[1]
  frame$Stage <- factor(train@misc$training_record_stage[observed],levels=c('PT','AT'))
  frame$Stage[is.na(frame$Stage)] <- 'PT'
  heterogeneous <- length(unique(train@misc$training_record_stage))>1L
  if(heterogeneous) frame <- frame[order(frame$Stage),,drop=FALSE]
  residual_formula <- if(heterogeneous) ~dsum(~units|Stage) else ~units
  inverse <- asreml::ainverse(ped)
  asreml::asreml.options(trace=FALSE,workspace=getOption('almond.asreml_workspace','16gb'))
  model <- asreml::asreml(fixed=Pheno~1+Year,random=~vm(Ind,inverse),
                         residual=residual_formula,na.action=asreml::na.method(y='include'),data=frame)
  for (attempt in seq_len(5)) {
    if (isTRUE(model$converge)) break
    model <- asreml::update.asreml(model)
  }
  if (!isTRUE(model$converge)) stop('Pedigree model failed to converge after five updates.')
  coefficients <- model$coef$random
  coefficient_names <- if(is.null(dim(coefficients))) names(coefficients) else rownames(coefficients)
  coefficient_ids <- sub('.*_','',coefficient_names)
  keep <- coefficient_ids %in% as.character(prediction@id)
  if (!any(keep)) stop('Cannot identify ASReml individual coefficients.')
  values <- setNames(as.numeric(coefficients)[keep],coefficient_ids[keep])
  attr(values,'training_diagnostics') <- list(solver='ASReml',heterogeneous=heterogeneous,
    record_counts=table(train@misc$training_record_stage),variance_components=summary(model)$varcomp,
    converged=model$converge,fixed_year=TRUE,cohort_keys=unique(train@misc$training_record_cohort))
  values
}

# Exact, dense pedigree BLUP with known simulated variance components, used
# only for small validation fixtures when an ASReml license is unavailable.
redesign_validation_blup <- function(state,train) {
  size <- SP$lastId
  if(size>3000) stop('Validation solver is limited to small fixtures; use ASReml for the experiment.')
  relationship <- matrix(0,size,size)
  for(i in seq_len(size)) {
    parents <- SP$pedigree[i,]
    known <- parents[parents>0]
    if(i>1 && length(known)) {
      relationship[i,seq_len(i-1)] <- colSums(relationship[known,seq_len(i-1),drop=FALSE])/2
      relationship[seq_len(i-1),i] <- relationship[i,seq_len(i-1)]
    }
    relationship[i,i] <- 1+if(length(known)==2) relationship[known[1],known[2]]/2 else 0
  }
  ids <- as.integer(train@id)
  response <- as.numeric(train@pheno)
  covariance <- state$reference_variance*relationship[ids,ids,drop=FALSE]
  residual_variances <- vapply(as.character(train@id),function(id) state$errors[[id]]$variance,numeric(1))
  diag(covariance) <- diag(covariance)+residual_variances
  design <- model.matrix(~factor(train@fixEff))
  solved <- solve(covariance,cbind(response,design))
  fixed <- solve(crossprod(design,solved[,-1,drop=FALSE]),crossprod(design,solved[,1]))
  weights <- solved[,1]-solved[,-1,drop=FALSE]%*%fixed
  values <- state$reference_variance*relationship[,ids,drop=FALSE]%*%weights
  setNames(as.numeric(values),as.character(seq_len(size)))
}

redesign_predict <- function(pop,method,model) {
  if (is.null(pop) || method=='Pheno') return(pop)
  if (method=='GS') return(setEBV(pop,model,value='bv',simParam=SP))
  values <- unname(model[as.character(pop@id)])
  if (anyNA(values)) stop('Missing pedigree predictions for selection candidates.')
  pop@ebv <- matrix(values,ncol=1)
  pop
}

redesign_choose <- function(pop,n,method,model=NULL) {
  pop <- redesign_predict(pop,method,model)
  selectInd(pop,nInd=n,use=if(method=='Pheno') 'pheno' else 'ebv',simParam=SP)
}

# Fixed annual replacement: new candidates are ranked separately from the
# current parents. Their merit does not determine whether the quota is filled.
redesign_renew <- function(parents,candidates,method,model,config) {
  quota <- as.integer(config$entrants)
  stopifnot(quota>=0L,quota<=config$parents,nInd(parents)==config$parents)
  eligible <- which(!candidates@id %in% parents@id)
  if(length(eligible)<quota) stop('Insufficient distinct new candidates for obligatory parental replacement.')
  newcomers <- if(quota>0L) redesign_choose(candidates[eligible],quota,method,model) else NULL
  retained <- if(quota<config$parents) {
    redesign_choose(parents,config$parents-quota,method,model)
  } else NULL
  chosen <- redesign_merge_latest(list(retained,newcomers))
  new <- sum(!chosen@id %in% parents@id)
  stopifnot(new==quota,nInd(chosen)==config$parents)
  list(parents=chosen,new=new)
}

redesign_advance <- function(state) {
 previous <- state$pop
 for(i in length(state$config$stages):2L) state$pop[[state$config$stages[i]]] <- previous[[state$config$stages[i-1L]]]
 state$pop$F1 <- NULL
}

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

redesign_year <- function(state,method,first_parental_use,rep,record=TRUE,renew=TRUE) {
 config <- state$config
 if(!config$design_version%in%c('PT_AT_ET_optiSel_single_genotype_Ne_v10','PT_AT_ET_optiSel_SI_groups_v11','PT_AT_ET_optiSel_adultcap32_v12')) stop('Historical burn-in is incompatible with the revised calendar.')
 first <- as.integer(first_parental_use)
 if(!first %in% c(9L,6L,4L,2L) || (method=='Pheno' && first<6L)) stop('Unsupported parental-use scheme.')
 state$training_stage <- config$training_stage[[as.character(first)]]
 state$training_cohort_stage <- config$training_cohort_stage[[as.character(first)]]
 target <- config$parent_stage[[as.character(first)]]
 redesign_advance(state)
 state$speed_cross_trace <- NULL;state$speed_ocs_solver <- NULL
 entrants <- 0L;allocation <- NULL
 if(renew && !is.null(state$pending)) {
  allocation <- state$pending
  if(allocation$use_year!=state$year) stop('Pending parental selection has the wrong crossing year.')
  state$parents <- allocation$parents;entrants <- allocation$new;state$pending <- NULL
 }
 candidates <- state$pop[[target]]
 start_model <- NULL
 if(method!='Pheno') {
  if(identical(state$model_method,method)) start_model <- state$prediction_pt_model
  if(is.null(start_model)) start_model <- redesign_fit_for_stage(state,method,'PT2')
 }
 accuracy <- NA_real_
 seedlings_before <- if(is.null(state$pop$Seedlings)) 0L else nInd(state$pop$Seedlings)
 if(seedlings_before && method!='Pheno') {
  predicted <- redesign_predict(state$pop$Seedlings,method,start_model)
  accuracy <- suppressWarnings(cor(as.numeric(gv(predicted)),as.numeric(ebv(predicted))))
 }
 if(isTRUE(config$early_selection) && seedlings_before) {
  if(method=='Pheno') stop('Phenotypic early selection is unavailable before yield records.')
  state$pop$Seedlings <- redesign_choose(state$pop$Seedlings,min(config$seedlings_keep,seedlings_before),method,start_model)
 }
 if(renew && first<6L && !is.null(candidates)) {
  # Own PT1 record is still missing; juvenile and year-4 selection are predictive.
  if(any(is.finite(as.numeric(candidates@pheno)))) stop('Pre-phenotyping candidates unexpectedly have yield records.')
  allocation <- redesign_allocate(state,candidates,method,start_model,first,state$year)
  state$parents <- allocation$parents;entrants <- allocation$new
 }
 ocs_solution <- if(is.null(allocation)) NULL else allocation$ocs
 options(almond.cross_state=state)
 on.exit(options(almond.cross_state=NULL),add=TRUE)
 if(is.null(ocs_solution)) {
  state$pop$F1 <- randCrossGamSI(list(SIPos=locuscompt_position_vector,random_failure=FALSE),
    state$parents,nCrosses=config$crosses,nProgeny=config$progeny)
 } else {
  state$speed_ocs_solver <- ocs_solution$solver_details
  state$pop$F1 <- randCrossGamSI(list(SIPos=locuscompt_position_vector,random_failure=FALSE),
    state$parents,nCrosses=NA_integer_,nProgeny=config$progeny,
    cross_plan=ocs_solution$cross_ids,family_sizes=ocs_solution$family_sizes)
  realized <- table(factor(c(state$pop$F1@mother,state$pop$F1@father),levels=state$parents@id))/(2*nInd(state$pop$F1))
  expected <- ocs_solution$contributions$contribution[match(names(realized),ocs_solution$contributions$id)]
  if(anyNA(expected)||max(abs(as.numeric(realized)-expected))>1e-10) stop('Realized OCS contributions differ from the validated plan.')
  state$ocs_trace[[as.character(state$year)]] <- ocs_solution
 }
 state$birth[as.character(state$pop$F1@id)] <- state$year
 if(nInd(state$pop$F1)!=config$crosses*config$progeny) stop('Incomplete compatible progeny cohort.')
 for(stage in setdiff(names(config$h2),'Parents')) state$pop[[stage]] <- redesign_measure(state$pop[[stage]],stage,state)
 redesign_measure_parents(state)
 redesign_archive_cohorts(state)
 pt_model <- redesign_fit_for_stage(state,method,'PT2')
 at_model <- redesign_fit_for_stage(state,method,'AT3')
 decision_models <- list(PT2=pt_model,AT3=at_model)
 end_model <- if(first==9L) at_model else pt_model
 state$prediction_pt_model <- pt_model;state$prediction_at_model <- at_model
 state$prediction_model <- pt_model;state$model_method <- method;state$model_training_stage <- 'PT2'
 next_selection <- NULL
 if(renew && first>=6L && !is.null(state$pop[[target]])) {
  if(method!='Pheno') redesign_assert_decision_training(state$pop[[target]],end_model,target)
  next_selection <- redesign_allocate(state,state$pop[[target]],method,end_model,first,state$year)
  state$pending <- next_selection
 }
 state$parent_accuracy <- if(!is.null(next_selection)) next_selection$accuracy else if(!is.null(allocation)) allocation$accuracy else NA_real_
 redesign_trial_cuts(state,method,decision_models)
 if(!record) return(NULL)
 scenario <- if(!is.null(config$scenario_id)) config$scenario_id else paste0(method,'_Y',first)
 row <- redesign_record(state,scenario,rep,accuracy,entrants,target)
 row$ocsMaxContribution <- if(is.na(config$ocs_alpha)) NA_real_ else config$ocs_max_contribution
 row$firstParentalUse <- first
 row$parentCandidateAge <- if(first>=6L) first-1L else first
 row$parentPhenotypingYears <- if(first==9L) 5L else if(first==6L) 2L else 0L
 row$parentSelectionYear <- if(is.null(allocation)) NA_integer_ else allocation$selected_year
 row$firstUseOfSelectedParents <- if(is.null(allocation)) NA_integer_ else allocation$use_year
 row$selectedForNextYear <- if(is.null(next_selection)) 0L else next_selection$new
 row$trialCutStage <- 'PT2';row$eliteCutStage <- 'AT3'
 row$predictionSolver <- if(method=='GS') if(first==9L) 'AlphaSimR::solveRRBLUP' else 'AlphaSimR::RRBLUP' else if(method=='Pedigree') 'ASReml' else NA_character_
 row$trainingPTRecords <- if(method=='Pheno') NA_integer_ else sum(redesign_training(state,candidates)@misc$training_record_stage=='PT')
 row$trainingATRecords <- if(method=='Pheno') NA_integer_ else if(first==9L) 75L else 0L
 row$trainingCohortAnchor <- if(method=='Pheno') NA_character_ else state$training_cohort_stage
 row$annualEnvironmentP <- state$environment_p[[as.character(state$year)]]
 row$trainingStage <- if(method=='Pheno') NA_character_ else state$training_stage
 row$trialTrainingAnchor <- if(method=='Pheno') NA_character_ else 'PT2'
 row$eliteTrainingAnchor <- if(method=='Pheno') NA_character_ else 'AT3'
 row$trialTrainingPTRecords <- if(method=='Pheno') NA_integer_ else 1500L
 row$eliteTrainingPTRecords <- if(method=='Pheno') NA_integer_ else 1500L
 row$trialTrainingATRecords <- if(method=='Pheno') NA_integer_ else 0L
 row$eliteTrainingATRecords <- if(method=='Pheno') NA_integer_ else if(first==9L) 75L else 0L
 row$trainingRecords <- if(method=='Pheno') NA_integer_ else if(first==9L) 1575L else 1500L
 row$trainingVariant <- if(method=='Pheno') NA_character_ else if(first==9L) 'Y9_PT1500_plus_AT75_calibrated' else 'decision_specific_PT1500'
 row$residualShape <- if(method=='Pheno') NA_character_ else if(first==9L) 'design_calibrated_PT_AT' else 'common_PT'
 row$ownPTRecordsValidated <- method!='Pheno'
 row$trainingCohorts <- if(method=='Pheno') NA_integer_ else min(state$config$training_cohorts,length(state$stage_history[[state$training_cohort_stage]]))
 row$trainingIndividuals <- if(method=='Pheno') NA_integer_ else nInd(redesign_training(state,candidates))
 row$earlySelection <- isTRUE(config$early_selection)
 row$progenyPerCross <- if(is.null(ocs_solution)) config$progeny else NA_integer_
 row$variableFamilySizes <- !is.null(ocs_solution)
 row$numberOfParentalPairs <- if(is.null(state$speed_cross_trace)) NA_integer_ else state$speed_cross_trace$distinct_pair_count
 row$maxContributingParents <- if(is.na(config$ocs_alpha)) config$parents else NA_integer_
 row$seedlingsBeforeSelection <- seedlings_before
 row$seedlingsAfterSelection <- nInd(state$pop$Seedlings)
 row$targetNewParents <- if(is.na(config$ocs_alpha)) config$entrants else NA_integer_
 row$targetRenewalPercent <- if(is.na(config$ocs_alpha)) 100*config$entrants/config$parents else NA_real_
 row$realizedRenewalPercent <- 100*entrants/nInd(state$parents)
 row$renewalRule <- if(is.na(config$ocs_alpha)) 'exact' else 'OCS_free'
 row$parentsUsed <- nInd(state$parents)
 row$meanParentalAge <- mean(state$year-state$birth[as.character(state$parents@id)])
 trace <- state$speed_cross_trace
 row$minimumPollenAge <- config$minimum_pollen_age;row$minimumSeedParentAge <- config$minimum_mother_age
 row$meanMotherAge <- mean(trace$mother_ages);row$meanFatherAge <- mean(trace$father_ages)
 row$weightedParentalAge <- (row$meanMotherAge+row$meanFatherAge)/2
 row$pollenFromAge2 <- mean(trace$father_ages==2L)
 row$pollenFromUnder4 <- mean(trace$father_ages<4L)
 row$ocsAlpha <- NA_real_;row$ocsTargetNe <- config$ocs_target_Ne;row$ocsImplementation <- if(is.na(config$ocs_alpha)) NA_character_ else 'optiSel_adultcap32_v12'
 row$ocsBound <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$bound
 row$ocsGenomicCoancestry <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$metrics$kinship
 row$ocsMateCoancestry <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$metrics$mate_kinship
 row$ocsExpectedBV <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$metrics$gain
 row$ocsEffectiveParents <- if(is.null(ocs_solution)) NA_real_ else 1/sum(ocs_solution$metrics$contribution^2)
 row$ocsBoundSatisfied <- if(is.null(ocs_solution)) NA else ocs_solution$metrics$kinship<=ocs_solution$bound+1e-8
 row$ocsCoancestryExcess <- if(is.null(ocs_solution)) NA_real_ else max(0,ocs_solution$metrics$kinship-ocs_solution$bound)
 row$ocsCoancestryExcessRelative <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$solver_details$optimization_audit$rounded_kinship_excess_relative
 row$ocsEndpointsCollapsed <- if(is.null(ocs_solution)) NA else abs(ocs_solution$maxK-ocs_solution$minK)<1e-10
 row$ocsFallbackToDiversity <- FALSE
 row$designVersion <- config$design_version
 row
}

# Refuse incomplete banks: 32 distinct diploids cannot be obtained from the
# historical files containing only 30 founders.
redesign_founder_bank <- function(bank_root) {
  paths <- list.files(bank_root,pattern='^run_',full.names=TRUE)
  paths <- paths[dir.exists(paths)]
  if(file.exists(file.path(bank_root,'manifest.json'))) paths <- c(bank_root,paths)
  if(length(paths)) paths <- paths[order(file.info(paths)$mtime,decreasing=TRUE)]
  for(path in paths) {
    manifest_file <- file.path(path,'manifest.json')
    if(!file.exists(manifest_file)) next
    manifest <- jsonlite::fromJSON(manifest_file)
    if(!identical(manifest$status,'complete')) next
    if(manifest$configuration$n_diploid_individuals!=150L ||
       manifest$configuration$mutation_rate!=1e-8 || nrow(manifest$chromosomes)!=8L ||
       !setequal(as.character(manifest$chromosomes$chromosome),as.character(1:8)) ||
       anyNA(manifest$chromosomes$exported_sites) ||
       any(manifest$chromosomes$exported_sites<=0)) next
    files <- c(paste0('Chr_',1:8,'_position.txt'),paste0('Chr_',1:8,'_gmatrix.txt'))
    if(all(file.exists(file.path(path,files)))) return(path)
  }
  stop('No complete bank of 150 diploid founders with eight chromosomes and mutation_rate=1e-8 is available yet. Wait for bank150 generation to finish.')
}

# A temporary smaller founder sample can establish a 32-parent program after
# pipeline fill. Added parents must already have eight phenotyping years.
redesign_fill_parent_pool <- function(state) {
 if(nInd(state$parents)!=state$config$parents) stop('The revised design requires exactly 32 sampled founders, not a legacy pool.')
 invisible(state)
}

redesign_setup <- function(root,config) {
  # Reuse validated founder import and S-locus assignment from the original program.
  source(file.path(root,'00_Burn_in','GlobalParameters.R'),local=.GlobalEnv)
  source(file.path(root,'compatible_crosses.R'),local=.GlobalEnv)
  source(file.path(root,'ocs_discrete.R'),local=.GlobalEnv)
  source(file.path(root,'manual_ocs_adapter.R'),local=.GlobalEnv)
  if(!identical(config$founder_source,'bank150')) stop('The revised design requires the bank150/60000 founder source.')
  bank <- redesign_founder_bank(config$founder_bank)
  founder_count <- config$parents
  assign('nParents',founder_count,envir=.GlobalEnv)
  options(almond.haplotype_dir=bank,almond.founder_count=founder_count)
  source(file.path(root,'CreateFounders_redesigned.R'),local=.GlobalEnv)
  stopifnot(nInd(Parents)==founder_count)
  state <- new.env(parent=emptyenv())
  state$config <- config
  state$founder_count <- founder_count
  state$founder_indices <- founder_indices
  state$ocs_frequency <- ocs_fixed_frequency(Parents)
  state$ocs_trace <- list()
  state$reference_variance <- initVarG
  state$reference_gxe_variance <- initVarGE
  state$annual_evaluations <- new.env(parent=emptyenv())
  environment_seed <- if(is.null(config$environment_seed)) config$seed+900000L else config$environment_seed
  state$environment_p <- redesign_environment_schedule(config,environment_seed)
  state$errors <- new.env(parent=emptyenv())
  state$parents <- Parents
  state$birth <- setNames(rep(-30L,nInd(Parents)),as.character(Parents@id))
  state$pop <- setNames(vector('list',length(config$stages)),config$stages)
  state$stage_history <- list()
  state$training_stage <- NULL
  state$pending <- NULL
  state$prediction_model <- NULL
  state$model_method <- NULL
  state$model_training_stage <- NULL
  state$year <- -config$fill_cycles
  state$parents <- redesign_measure(state$parents,'Parents',state)
  for(i in seq_len(config$fill_cycles)) {
    state$year <- i-config$fill_cycles
    redesign_year(state,'Pheno',9L,0L,record=FALSE,renew=FALSE)
  }
  redesign_fill_parent_pool(state)
  state
}

redesign_bind <- function(frames) {
 frames<-Filter(function(x) !is.null(x)&&nrow(x)>0L,frames)
 if(!length(frames)) return(NULL)
 columns<-unique(unlist(lapply(frames,names)))
 do.call(rbind,lapply(frames,function(x) {
   for(column in setdiff(columns,names(x))) x[[column]]<-NA
   x[columns]
 }))
}

redesign_attempt_scenario <- function(action) {
 tryCatch(list(ok=TRUE,result=action()),error=function(e)
   list(ok=FALSE,error=conditionMessage(e)))
}

redesign_all_scenarios <- function(root=NULL) redesign_scenarios(root)

# Cohorts are archived before either fixed trial cut. Each ID contributes only
# its latest aggregate at the selected training stage; parental updates are excluded.
redesign_trial_cuts <- function(state,method,models) {
 for(entry in list(c('PT2',state$config$act),c('AT3',state$config$ect))) {
  stage <- entry[1];n <- as.integer(entry[2]);pop <- state$pop[[stage]]
  if(!is.null(pop)) {
   model <- if(method=='Pheno') NULL else models[[stage]]
   if(method!='Pheno') redesign_assert_decision_training(pop,model,stage)
   state$pop[[stage]] <- redesign_choose(pop,min(n,nInd(pop)),method,model)
  }
 }
 invisible(state)
}

redesign_worker <- function(job) {
 .libPaths(job$libraries)
 options(almond.input_dir=job$root,almond.nThreads=job$config$threads)
 source(file.path(job$root,'redesigned_schemes.R'),local=.GlobalEnv)
 assign('.Random.seed',job$seed,envir=.GlobalEnv)
 job$config$environment_seed <- job$config$seed+900000L+as.integer(job$rep)
 dir.create(job$directory,recursive=TRUE,showWarnings=FALSE)
 timing_file<-file.path(job$directory,paste0('timings_',job$phase,'.tsv'))
 options(almond.timing_file=timing_file,almond.timing_rep=job$rep,
   almond.timing_scenario='initialization',almond.timing_year=0L,
   almond.diagnostic_file=file.path(job$directory,paste0('diagnostics_',job$phase,'.log')))
 install_timers<-function() {
   for(name in c('redesign_setup','redesign_advance','redesign_measure','redesign_fit',
     'redesign_predict','redesign_choose','redesign_renew','redesign_ocs',
     'randCrossGamSI','redesign_record','redesign_year','redesign_allocate','redesign_measure_parents')) {
     if(exists(name,envir=.GlobalEnv,inherits=FALSE))
       assign(name,redesign_timer(get(name,envir=.GlobalEnv),name),envir=.GlobalEnv)
   }
 }
 status_file<-file.path(job$directory,paste0('progress_',job$phase,'.json'))
 progress<-function(status,scenario,year=0L,error=NULL) {
   options(almond.timing_scenario=scenario,almond.timing_year=year)
   jsonlite::write_json(list(status=status,rep=job$rep,phase=job$phase,
     scenario=scenario,year=year,memory=redesign_memory(),
     updated=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),error=error),
     status_file,pretty=TRUE,auto_unbox=TRUE)
 }
 log<-file(file.path(job$directory,paste0('execution_',job$phase,'.log')),'wt')
 sink(log);on.exit({sink();close(log)},add=TRUE)
 base_file<-file.path(job$directory,'burnin_redesigned.RData')
 initialized<-redesign_attempt_scenario(function() {
   progress('running','initialization')
   if(job$phase=='parallel') {
     install_timers()
     state<-redesign_setup(job$root,job$config)
     burnin<-vector('list',job$config$burnin)
     for(year in seq_len(job$config$burnin)) {
       state$year<-year;progress('running','Burn_in',year)
       cat('Burn-in year',year,'\n')
       burnin[[year]]<-redesign_year(state,'Pheno',9L,job$rep)
       burnin[[year]]$scenario<-'Burn_in'
     }
     save(state,burnin,SP,locuscompt_position_vector,self_compatible_allele,file=base_file)
     saveRDS(redesign_bind(burnin),file.path(job$directory,'burnin_results.rds'))
     saveRDS(state$environment_p,file.path(job$directory,'annual_environment.rds'))
   } else if(!file.exists(base_file)) stop('This replicate has no successfully completed burn-in.')
   TRUE
 })
 if(!initialized$ok) {
   progress('failed','Burn_in',getOption('almond.timing_year'),initialized$error)
   writeLines(initialized$error,file.path(job$directory,paste0('error_',job$phase,'.txt')))
   status<-list(ok=FALSE,rep=job$rep,error=initialized$error,fatal_initialization=TRUE)
   saveRDS(status,file.path(job$directory,paste0('worker_status_',job$phase,'.rds')))
   return(status)
 }
 outputs<-list();statuses<-list()
 selected<-job$scenarios[job$scenarios$method%in%job$methods,,drop=FALSE]
 for(index in seq_len(nrow(selected))) {
   scenario<-selected$scenario_id[index]
   records<-list();completed_year<-job$config$burnin
   started<-Sys.time()
   outcome<-redesign_attempt_scenario(function() {
     # Reinstall clean functions and load a fresh paired starting population.
     # A failed scenario must not contaminate the next scenario's state or RNG.
     source(file.path(job$root,'redesigned_schemes.R'),local=.GlobalEnv)
     source(file.path(job$root,'ocs_discrete.R'),local=.GlobalEnv)
     source(file.path(job$root,'manual_ocs_adapter.R'),local=.GlobalEnv)
     source(file.path(job$root,'compatible_crosses.R'),local=.GlobalEnv)
     load(base_file,envir=.GlobalEnv)
     state<-get('state',.GlobalEnv)
     if(!identical(state$config$design_version,job$config$design_version) && !(state$config$design_version %in% c('PT_AT_ET_annual_matched_PT_v4','PT_AT_ET_decision_specific_PT_v5','PT_AT_ET_optiSel_profiles_v7','PT_AT_ET_optiSel_cap10_v8','PT_AT_ET_optiSel_native_offspring_v9') && job$config$design_version%in%c('PT_AT_ET_optiSel_single_genotype_Ne_v10','PT_AT_ET_optiSel_SI_groups_v11','PT_AT_ET_optiSel_adultcap32_v12'))) stop('Historical burn-in cannot be reused with the revised calendar.')
     if(is.null(state$stage_history)) stop('Burn-in lacks stage-specific training archives.')
     SP$nThreads<-job$config$threads
     assign('.Random.seed',parallel::nextRNGStream(job$seed),envir=.GlobalEnv)
     method<-selected$method[index];years<-selected$first_parental_use[index]
     state$config$ocs_alpha<-selected$ocs_alpha[index]
     state$config$ocs_target_Ne<-selected$ocs_target_Ne[index]
     state$config$early_selection<-selected$early_selection[index]
     state$config$progeny<-selected$progeny_per_cross[index]
     state$config$seedlings_keep<-selected$seedlings_keep[index]
     state$config$entrants<-if(is.na(selected$renewal_percent[index]))
       selected$new_parents_per_year[index] else as.integer(state$config$parents*selected$renewal_percent[index]/100)
     state$config$scenario_id<-scenario
     state$config$speed_breeding<-isTRUE(selected$speed_breeding[index])
     state$config$minimum_pollen_age<-if(state$config$speed_breeding) 2L else 4L
     state$config$minimum_mother_age<-4L
     state$config$first_parental_use<-years
     state$config$ocs_backend<-'optiSel_adultcap32_v12'
     state$config$ocs_max_contribution<-job$config$ocs_max_contribution
     state$config$ocs_family_sizes<-job$config$ocs_family_sizes
     state$config$ocs_max_parents<-NA_integer_
     if(!identical(state$config$design_version,job$config$design_version) && !(state$config$design_version %in% c('PT_AT_ET_annual_matched_PT_v4','PT_AT_ET_decision_specific_PT_v5','PT_AT_ET_optiSel_profiles_v7','PT_AT_ET_optiSel_cap10_v8','PT_AT_ET_optiSel_native_offspring_v9') && job$config$design_version%in%c('PT_AT_ET_optiSel_single_genotype_Ne_v10','PT_AT_ET_optiSel_SI_groups_v11','PT_AT_ET_optiSel_adultcap32_v12'))) stop('Burn-in design version differs. Generate a new burn-in for matched-cohort training.')
     state$config$training_stage<-job$config$training_stage
     state$config$training_cohort_stage<-job$config$training_cohort_stage
     state$config$gs_solver<-job$config$gs_solver
     state$config$design_version<-job$config$design_version
     state$prediction_pt_model<-NULL;state$prediction_at_model<-NULL
     options(almond.decision_audit_dir=file.path(job$directory,'training_audit',scenario))
     state$prediction_model<-NULL;state$model_method<-NULL
     # Re-select at the end of the new scenario's year; do not import the
     # pending baseline choice or another scenario's genomic model.
     state$pending<-NULL
     if(method=='Pedigree') {
       if(!requireNamespace('asreml',quietly=TRUE)) stop('ASReml is unavailable.')
       asreml::asreml.options(workspace=job$config$asreml_workspace,trace=FALSE)
     }
     install_timers()
     cat('Scenario',scenario,'\n')
     for(k in seq_len(job$config$future)) {
       state$year<-job$config$burnin+k
       progress('running',scenario,state$year)
       records[[k]]<<-redesign_year(state,method,years,job$rep)
       records[[k]]$scenario<<-scenario
       completed_year<<-state$year
       partial<-redesign_bind(records)
       saveRDS(partial,file.path(job$directory,paste0('results_partial_',scenario,'.rds')))
       if(!is.na(state$config$ocs_alpha)||state$config$speed_breeding) {
         saveRDS(list(solver=state$speed_ocs_solver,crosses=state$speed_cross_trace),
           file.path(job$directory,sprintf('allocation_%s_year_%03d.rds',scenario,state$year)))
       }
       cat('Completed future year',k,'of',job$config$future,'\n')
     }
     result<-redesign_bind(records)
     stopifnot(nrow(result)==job$config$future,!anyDuplicated(result[c('rep','scenario','year')]))
     saveRDS(data.frame(id=seq_len(SP$lastId),dam=SP$pedigree[,1],sire=SP$pedigree[,2]),
       file.path(job$directory,paste0('pedigree_',scenario,'.rds')))
     saveRDS(result,file.path(job$directory,paste0('results_',scenario,'.rds')))
     result
   })
   statuses[[scenario]]<-data.frame(rep=job$rep,phase=job$phase,scenario=scenario,
     status=if(outcome$ok) 'complete' else 'failed',completed_years=completed_year-job$config$burnin,
     error=if(outcome$ok) NA_character_ else outcome$error,
     elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),
     stringsAsFactors=FALSE)
   saveRDS(redesign_bind(statuses),file.path(job$directory,paste0('scenario_status_',job$phase,'.rds')))
   if(outcome$ok) {
     outputs[[scenario]]<-outcome$result
     progress('scenario_complete',scenario,completed_year)
   } else {
     progress('scenario_failed',scenario,completed_year+1L,outcome$error)
     writeLines(outcome$error,file.path(job$directory,paste0('error_',scenario,'.txt')))
     cat('SCENARIO FAILED; continuing other schemes:',scenario,outcome$error,'\n')
   }
   saveRDS(redesign_bind(outputs),file.path(job$directory,paste0('results_',job$phase,'.rds')))
 }
 manifest<-redesign_bind(statuses)
 failed<-if(is.null(manifest)) FALSE else any(manifest$status=='failed')
 result<-list(ok=TRUE,rep=job$rep,result=redesign_bind(outputs),
   scenario_status=manifest,has_failures=failed)
 progress(if(failed) 'finished_with_failures' else 'complete','all')
 saveRDS(result[names(result)!='result'],file.path(job$directory,paste0('worker_status_',job$phase,'.rds')))
 result
}

run_redesigned_schemes <- function(root,config=redesign_config(),scenarios=redesign_scenarios()) {
 experiment_start<-proc.time()[['elapsed']]
 original_root<-normalizePath(root,winslash='/',mustWork=TRUE)
 if(!identical(config$founder_source,'bank150')) stop('The revised design requires bank150.')
 redesign_founder_bank(config$founder_bank)
 config$data_dir<-redesign_data_dir(original_root,config)
 if(is.null(config$asreml_workspace)) config$asreml_workspace<-'16gb'
 parent<-file.path(config$data_dir,'manuscript_runs')
 dir.create(parent,recursive=TRUE,showWarnings=FALSE)
 run<-tempfile(paste0('PT_AT_ET_',config$reps,'reps_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),tmpdir=parent)
 dir.create(run,recursive=TRUE)
 # Freeze source so later manuscript/code edits cannot change a running run.
 root<-file.path(run,'code');dir.create(root)
 dir.create(file.path(root,'00_Burn_in'))
 files<-list.files(original_root,pattern='[.]R$',full.names=TRUE)
 stopifnot(all(file.copy(files,root)))
 stopifnot(file.copy(file.path(original_root,'00_Burn_in','GlobalParameters.R'),file.path(root,'00_Burn_in')))
 writeLines(run,file.path(config$data_dir,'latest_run.txt'))
 RNGkind("L'Ecuyer-CMRG");set.seed(config$seed);stream<-.Random.seed
 jobs<-vector('list',config$reps)
 for(rep in seq_len(config$reps)) {
   jobs[[rep]]<-list(rep=rep,root=root,config=config,scenarios=scenarios,
     methods=c('Pheno','GS'),phase='parallel',seed=stream,libraries=.libPaths(),
     directory=file.path(run,sprintf('rep_%03d',rep)))
   dir.create(jobs[[rep]]$directory)
   saveRDS(jobs[[rep]],file.path(jobs[[rep]]$directory,'job_parallel.rds'))
   stream<-parallel::nextRNGStream(stream)
 }
 saveRDS(list(config=config,scenarios=scenarios,source_root=original_root),file.path(run,'configuration.rds'))
 write.csv(scenarios,file.path(run,'scenarios.csv'),row.names=FALSE)
 jsonlite::write_json(list(status='running',phase='parallel',replicates=config$reps,
   schemes=nrow(scenarios),started=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')),
   file.path(run,'experiment_progress.json'),auto_unbox=TRUE,pretty=TRUE)
 cluster<-parallel::makeCluster(min(config$workers,config$reps),outfile=file.path(run,'cluster.log'))
 on.exit(parallel::stopCluster(cluster),add=TRUE)
 results<-parallel::parLapplyLB(cluster,jobs,function(job) {
   .libPaths(job$libraries);source(file.path(job$root,'redesigned_schemes.R'));redesign_worker(job)
 })
 # Scheme failures do not stop unrelated schemes, replicates or pedigree.
 jsonlite::write_json(list(status='running',phase='pedigree',replicates=config$reps,
   schemes=nrow(scenarios)),file.path(run,'experiment_progress.json'),auto_unbox=TRUE,pretty=TRUE)
 pedigree<-vector('list',config$reps)
 for(rep in seq_len(config$reps)) {
   job<-jobs[[rep]];job$methods<-'Pedigree';job$phase<-'pedigree'
   pedigree[[rep]]<-parallel::parLapply(cluster[1],list(job),function(job) {
     .libPaths(job$libraries);source(file.path(job$root,'redesigned_schemes.R'));redesign_worker(job)
   })[[1]]
 }
 outcomes<-c(results,pedigree)
 combined<-redesign_bind(lapply(outcomes,function(x) x$result))
 statuses<-redesign_bind(lapply(outcomes,function(x) x$scenario_status))
 fatal<-Filter(function(x) !isTRUE(x$ok),outcomes)
 complete<-!length(fatal)&&(nrow(scenarios)==0L || (!is.null(statuses)&&
   all(statuses$status=='complete')&&nrow(statuses)==config$reps*nrow(scenarios)))
 expected_rows<-config$reps*nrow(scenarios)*config$future
 if(complete && expected_rows>0L) stopifnot(nrow(combined)==expected_rows,!anyDuplicated(combined[c('rep','scenario','year')]))
 saveRDS(combined,file.path(run,'results_all_schemes.rds'))
 saveRDS(statuses,file.path(run,'scenario_status_all.rds'))
 saveRDS(list(complete=complete,finished=TRUE,replicates=config$reps,schemes=nrow(scenarios),
   expected_rows=expected_rows,verified_rows=if(is.null(combined)) 0L else nrow(combined),
   failed_scenarios=if(is.null(statuses)) NULL else statuses[statuses$status=='failed',],
   initialization_failures=lapply(fatal,function(x) x[c('rep','error')])),file.path(run,'status.rds'))
 saveRDS(list(elapsed_seconds=proc.time()[['elapsed']]-experiment_start,completed=Sys.time()),
   file.path(run,'total_runtime.rds'))
 jsonlite::write_json(list(status=if(complete) 'complete' else 'finished_with_failures',
   replicates=config$reps,schemes=nrow(scenarios),expected_rows=expected_rows,
   verified_rows=if(is.null(combined)) 0L else nrow(combined)),
   file.path(run,'experiment_progress.json'),auto_unbox=TRUE,pretty=TRUE)
 message('Experiment finished: ',run,'; all scenarios complete: ',complete)
 invisible(list(run_dir=run,data=combined,statuses=statuses))
}


redesign_data_dir <- function(root,config=list()) {
  configured <- config$data_dir
  if(is.null(configured) || !nzchar(configured)) configured <- Sys.getenv('ALMONDSIM_DATA_DIR',unset='')
  if(!nzchar(configured)) configured <- file.path(dirname(dirname(normalizePath(root,winslash='/',mustWork=TRUE))),'almondsim_data')
  if(!grepl('^([A-Za-z]:[/\\\\]|/)',configured)) configured <- file.path(getwd(),configured)
  resolved <- normalizePath(configured,winslash='/',mustWork=FALSE)
  source_root <- normalizePath(root,winslash='/',mustWork=TRUE)
  repository <- if(basename(source_root)=='02_ClonalBreeding') dirname(source_root) else source_root
  if(tolower(resolved)==tolower(repository)||startsWith(tolower(resolved),paste0(tolower(repository),'/'))) {
    stop('Simulation data must be stored outside the source repository.')
  }
  resolved
}




# Inclusive elapsed times: nested function durations must not be added together.
redesign_timer <- function(original,label) {
  force(original);force(label)
  function(...) {
    started <- Sys.time(); cpu <- proc.time(); success <- FALSE
    on.exit({
      destination <- getOption('almond.timing_file')
      if(!is.null(destination)) {
        duration <- proc.time()-cpu
        row <- data.frame(rep=getOption('almond.timing_rep'),
          scenario=getOption('almond.timing_scenario'),year=getOption('almond.timing_year'),
          task=label,started=format(started,'%Y-%m-%dT%H:%M:%S%z'),
          elapsed_seconds=unname(duration['elapsed']),
          cpu_seconds=unname(duration['user.self']+duration['sys.self']),success=success)
        memory <- redesign_memory()
        row$working_set_gb <- memory$working_set_gb
        row$private_gb <- memory$private_gb
        row$peak_private_gb <- memory$peak_private_gb
        tryCatch(write.table(row,destination,sep='\t',row.names=FALSE,
          col.names=!file.exists(destination),append=file.exists(destination),quote=FALSE),
          error=function(e) cat('Timing log could not be written:',conditionMessage(e),'\n'))
      }
    },add=TRUE)
    result <- withCallingHandlers(original(...),
      warning=function(w) redesign_diagnostic(label,'warning',w),
      error=function(e) redesign_diagnostic(label,'error',e))
    success <- TRUE;result
  }
}

redesign_check_asreml <- function() {
  # A separate preflight process releases its license before pedigree starts.
  script <- tempfile(fileext='.R');on.exit(unlink(script),add=TRUE)
  libraries <- paste(capture.output(dput(.libPaths())),collapse='')
  writeLines(c(paste0('.libPaths(',libraries,')'),
    "library(asreml)",
    "d <- data.frame(y=c(1,2,3,5,4,6),g=factor(rep(1:3,2)))",
    "invisible(asreml::asreml(y~1,random=~g,data=d,trace=FALSE))"),script)
  code <- system2(file.path(R.home('bin'),'Rscript.exe'),shQuote(script))
  if(code!=0L) stop('ASReml preflight failed before starting the experiment.')
  invisible(TRUE)
}



redesign_memory <- function() {
  result <- list(working_set_gb=NA_real_,private_gb=NA_real_,peak_private_gb=NA_real_)
  if(requireNamespace('ps',quietly=TRUE)) {
    memory <- tryCatch(ps::ps_memory_info(),error=function(e) NULL)
    if(!is.null(memory)) {
      result$working_set_gb <- unname(memory['rss'])/1024^3
      result$private_gb <- unname(memory['mem_private'])/1024^3
      result$peak_private_gb <- unname(memory['peak_pagefile'])/1024^3
    }
  }
  result
}

redesign_diagnostic <- function(task,kind,condition) {
  destination <- getOption('almond.diagnostic_file')
  if(is.null(destination)) return(invisible(NULL))
  report <- c(paste(format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),kind,task),
    conditionMessage(condition),paste('Call:',paste(deparse(conditionCall(condition)),collapse=' ')),
    paste('Memory:',paste(names(redesign_memory()),unlist(redesign_memory()),collapse='; ')),
    'Call stack:',vapply(sys.calls(),function(x) paste(deparse(x),collapse=' '),character(1)))
  tryCatch(cat(paste(report,collapse='\n'),'\n',file=destination,append=TRUE),
    error=function(e) cat('Diagnostic:',conditionMessage(condition),'\n'))
  invisible(NULL)
}


redesign_measure_parents <- function(state) {
 age <- state$year-state$birth[as.character(state$parents@id)]
 if(anyNA(age)) stop('Unknown parental birth year.')
 parts <- list();juvenile <- which(age<4L)
 if(length(juvenile)) {
  pop <- state$parents[juvenile];pop@pheno[,] <- NA_real_;parts[[1L]] <- pop
 }
 stages <- setdiff(names(state$config$h2),'Parents')
 counts <- pmin(length(stages),pmax(1L,age-3L))
 for(count in sort(unique(counts[age>=4L]))) {
  pop <- state$parents[which(age>=4L & counts==count)]
  parts[[length(parts)+1L]] <- redesign_measure(pop,stages[count],state)
 }
 state$parents <- mergePops(parts)
 invisible(state)
}


redesign_allocate <- function(state,candidates,method,model,first,selection_year) {
 old <- state$parents
 predicted <- redesign_predict(candidates,method,model)
 score <- if(method=='Pheno') predicted@pheno else predicted@ebv
 accuracy <- suppressWarnings(cor(as.numeric(gv(predicted)),as.numeric(score)))
 if(!is.na(state$config$ocs_alpha)) {
  solution <- redesign_ocs(state,candidates,model)
  updated <- list(parents=solution$parents,new=solution$new,ocs=solution)
 } else {
  solution <- NULL
  updated <- redesign_renew(old,candidates,method,model,state$config)
  updated$ocs <- NULL
 }
 updated$selected_year <- selection_year
 updated$use_year <- if(first>=6L) selection_year+1L else selection_year
 updated$candidate_age <- if(first>=6L) first-1L else first
 updated$candidate_ids <- as.character(candidates@id)
 updated$new_ids <- setdiff(as.character(updated$parents@id),as.character(old@id))
 updated$accuracy <- accuracy
 updated
}


# Y9 repeated-record production extension. Other schemes keep their PT-only
# model implementations. Both decision models always include candidates' own PT.
redesign_training_PT_original <- redesign_training
redesign_fit_PT_original <- redesign_fit
redesign_fit_for_stage_PT_original <- redesign_fit_for_stage
redesign_predict_PT_original <- redesign_predict
redesign_uses_repeated_Y9 <- function(state) isTRUE(state$config$first_parental_use==9L)

redesign_training <- function(state,candidates,anchor=NULL) {
  stage <- if(is.null(anchor)) state$training_cohort_stage else anchor
  train <- redesign_training_PT_original(state,candidates,anchor=stage)
  extra <- data.frame(id=character(),pheno=numeric(),record_year=integer(),stage=character(),cohort=character())
  if(redesign_uses_repeated_Y9(state) && stage=='AT3') {
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

redesign_repeated_records <- function(train) {
  pt <- data.frame(id=as.character(train@id),pheno=as.numeric(train@pheno),record_year=as.integer(train@fixEff),stage='PT',cohort=as.character(train@misc$training_record_cohort))
  records <- rbind(pt,train@misc$extra_AT_records)
  stopifnot(sum(records$stage=='PT')==1500L,length(unique(records$id))==1500L,!anyDuplicated(records[c('id','stage')]))
  if(nrow(train@misc$extra_AT_records))stopifnot(nrow(records)==1575L,sum(table(records$id)==2L)==75L,sum(records$stage=='AT')==75L)
  records
}

redesign_repeated_frame <- function(records,prediction) {
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

redesign_fit_repeated_pedigree <- function(state,method,candidates) {
  if(method=='Pheno')return(NULL)
  if(method!='Pedigree')stop('Repeated pedigree fitter requires method Pedigree.')
  if(state$training_cohort_stage=='PT2')return(redesign_fit_PT_original(state,method,candidates))
  train <- redesign_training(state,candidates)
  records <- redesign_repeated_records(train)
  prediction <- redesign_merge_latest(c(list(train,state$parents),state$pop))
  frame <- redesign_repeated_frame(records,prediction)
  ped <- data.frame(id=as.character(seq_len(SP$lastId)),dam=as.character(SP$pedigree[,1]),sire=as.character(SP$pedigree[,2]))
  ped$dam[ped$dam=='0']<-NA;ped$sire[ped$sire=='0']<-NA
  inverse <- asreml::ainverse(ped)
  asreml::asreml.options(trace=FALSE,workspace=getOption('almond.asreml_workspace','16gb'))
  model <- redesign_asreml_calibrated(state,frame,inverse)
  stable <- function(fit) {
    if(!isTRUE(fit$converge))return(FALSE)
    vc<-summary(fit)$varcomp
    change<-vc[['%ch']][vc$bound!='F']
    all(is.finite(change))&&all(abs(change)<=1)
  }
  for(attempt in seq_len(5)) {
    if(stable(model))break
    model<-asreml::update.asreml(model)
  }
  if(!stable(model))stop('Repeated-record model did not converge with stable variance/covariance components after five updates.')
  coefficients<-model$coef$random
  coefficient_names<-if(is.null(dim(coefficients)))names(coefficients)else rownames(coefficients)
  coefficient_ids<-sub('.*_','',coefficient_names)
  keep<-grepl('vm(',coefficient_names,fixed=TRUE)&coefficient_ids%in%as.character(prediction@id)
  if(!any(keep)||anyDuplicated(coefficient_ids[keep]))stop('Cannot uniquely identify pedigree individual effects.')
  values<-setNames(as.numeric(coefficients)[keep],coefficient_ids[keep])
  if(!all(as.character(prediction@id)%in%names(values)))stop('Missing pedigree effects for prediction IDs.')
  components<-summary(model)$varcomp
  component<-function(label)components$component[grepl(label,rownames(components),fixed=TRUE)]
  rPT<-component('Stage_PT:PT');rAT<-component('Stage_AT:AT');rCov<-component('Stage_AT:PT')
  if(length(rPT)!=1L||length(rAT)!=1L||length(rCov)!=1L)stop('Cannot identify the repeated residual covariance components.')
  R<-matrix(c(rPT,rCov,rCov,rAT),2,2,dimnames=list(c('PT','AT'),c('PT','AT')))
  if(any(!is.finite(R))||min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)<=0)stop('Repeated residual matrix is not positive definite.')
  calibration <- redesign_calibrated_residual(state)
  if(max(abs(R/R[1,1]-calibration$relative_matrix))>1e-7)stop('Pedigree residual ratios differ from the prespecified GS calibration.')
  attr(values,'training_diagnostics')<-list(solver='ASReml',heterogeneous=TRUE,repeated_residual=TRUE,
    residual_structure='id(Ind):us(Stage) constrained shape',residual_shape_source='design_calibrated',calibration=calibration,observed_records=nrow(records),unique_individuals=length(unique(records$id)),
    paired_individuals=sum(table(records$id)==2L),record_counts=table(records$stage),
    variance_components=components,residual_matrix=R,residual_correlation=rCov/sqrt(rPT*rAT),
    stable_components=stable(model),converged=model$converge,fixed_year=TRUE,cohort_keys=unique(records$cohort))
  values
}

redesign_fit_for_stage <- function(state,method,anchor) {
  if(method=='Pheno')return(NULL)
  if(!redesign_uses_repeated_Y9(state))return(redesign_fit_for_stage_PT_original(state,method,anchor))
  old_anchor<-state$training_cohort_stage;old_stage<-state$training_stage
  on.exit({state$training_cohort_stage<-old_anchor;state$training_stage<-old_stage},add=TRUE)
  state$training_cohort_stage<-anchor;state$training_stage<-'PT2'
  model<-redesign_fit(state,method,state$pop[[anchor]])
  train<-redesign_training(state,NULL,anchor=anchor);records<-redesign_repeated_records(train)
  meta<-list(anchor=anchor,ids=as.character(train@id),cohort_keys=unique(records$cohort),
    record_ids=records$id,record_years=records$record_year,record_stages=records$stage,records=records,
    PT_records=sum(records$stage=='PT'),AT_records=sum(records$stage=='AT'),total_records=nrow(records),
    diagnostics=attr(model,'training_diagnostics'))
  attr(model,'decision_training')<-meta
  audit<-getOption('almond.decision_audit_dir')
  if(!is.null(audit)) {
    dir.create(audit,recursive=TRUE,showWarnings=FALSE)
    saveRDS(meta,file.path(audit,sprintf('fit_year_%03d_%s_last_%s.rds',state$year,anchor,tail(meta$cohort_keys,1))))
  }
  model
}

redesign_calibrated_residual <- function(state) {
  A <- state$reference_variance
  GE <- state$reference_gxe_variance
  stages <- c('PT2','AT1','AT2','AT3')
  durations <- 2:5
  total <- A*(1/state$config$h2[stages]-1)
  pure <- as.numeric(total-GE/durations)
  stopifnot(all(is.finite(pure)),all(pure>0))
  # Match redesign_error_step at EACH intermediate annual aggregate. Marginal
  # pure error increases from AT1 to AT3 because the GE/n budget decreases.
  pure_cov <- pure[1]
  for(j in 2:4) pure_cov <- pure_cov*min(pure[j-1],pure[j])/pure[j-1]
  # Expected annual GE covariance: two common years out of 2 and 5.
  # This is a design-level marginal working covariance, not a fit to true gv
  # or genotype-specific GxE slopes / realised annual effects.
  total_cov <- pure_cov+GE/5
  R <- matrix(c(total[1],total_cov,total_cov,total[4]),2,2,
              dimnames=list(c('PT','AT'),c('PT','AT')))
  stopifnot(all(is.finite(R)),min(eigen(R,symmetric=TRUE,only.values=TRUE)$values)>0)
  list(absolute_working_matrix=R,relative_matrix=R/R[1,1],
       residual_correlation=R[1,2]/sqrt(R[1,1]*R[2,2]),
       intermediate_pure_error=pure,intermediate_stages=stages,
       pure_error_covariance=pure_cov,expected_GE_covariance=GE/5,
       source='prespecified simulation design; no true individual genetic values')
}

redesign_native_repeated_gs <- function(state,train) {
  records <- redesign_repeated_records(train)
  is_at <- records$stage=='AT'
  prediction <- redesign_merge_latest(c(list(train,state$parents),state$pop))
  genotype_train <- AlphaSimR::pullSnpGeno(train,simParam=SP)-1
  positions <- match(records$id,as.character(train@id))
  stopifnot(!anyNA(positions))
  M <- genotype_train[positions,,drop=FALSE]
  X <- model.matrix(~factor(record_year),data=records)
  y <- matrix(records$pheno,ncol=1)
  stopifnot(qr(X)$rank==ncol(X),all(is.finite(y)),all(is.finite(M)))
  calibration <- redesign_calibrated_residual(state)
  R <- calibration$relative_matrix
  # PT rows have unit relative variance. Whiten each observed AT conditional
  # on its original PT: eAT - cov/PTvar * ePT, then divide by conditional SD.
  if(any(is_at)) {
    at <- which(is_at);pt <- match(records$id[at],records$id[!is_at])
    stopifnot(all(!is_at[pt]),all(records$id[pt]==records$id[at]))
    b <- R['AT','PT']/R['PT','PT']
    conditional_sd <- sqrt(R['AT','AT']-R['AT','PT']^2/R['PT','PT'])
    y[at,] <- (y[at,,drop=FALSE]-b*y[pt,,drop=FALSE])/conditional_sd
    X[at,] <- (X[at,,drop=FALSE]-b*X[pt,,drop=FALSE])/conditional_sd
    M[at,] <- (M[at,,drop=FALSE]-b*M[pt,,drop=FALSE])/conditional_sd
  }
  fit <- AlphaSimR::solveRRBLUP(y,X,M)
  stopifnot(all(is.finite(fit$u)),is.finite(fit$Vu),fit$Vu>0,
            is.finite(fit$Ve),fit$Ve>0)
  genotype_prediction <- AlphaSimR::pullSnpGeno(prediction,simParam=SP)-1
  # Same marker effects on all candidate genotypes, with a constant centring
  # term that does not change ranks. One genetic effect per unique genotype.
  value <- drop(genotype_prediction%*%fit$u)-sum(colMeans(genotype_train)*drop(fit$u))
  values <- setNames(value,as.character(prediction@id))
  stopifnot(!anyDuplicated(names(values)),all(is.finite(values)))
  fitted_y <- X%*%fit$beta+M%*%fit$u
  residual <- y-fitted_y
  # Independent normal equation checks, including the residual scale.
  normal_fixed <- max(abs(crossprod(X,residual)))
  normal_markers <- max(abs(crossprod(M,residual)-(fit$Ve/fit$Vu)*fit$u))
  relative_error <- max(normal_fixed,normal_markers)/max(1,max(abs(crossprod(M,y))))
  if(!is.finite(relative_error)||relative_error>1e-6)stop('Weighted GS normal-equation validation failed.')
  diagnostic <- list(solver='AlphaSimR::solveRRBLUP',version=as.character(packageVersion('AlphaSimR')),
    repeated_residual=any(is_at),heterogeneous=any(is_at),fixed_year=TRUE,
    observed_records=nrow(records),unique_individuals=length(unique(records$id)),
    paired_individuals=sum(table(records$id)==2L),record_counts=table(records$stage),
    residual_shape_source='design_calibrated',calibration=calibration,
    residual_matrix=if(any(is_at))R*as.numeric(fit$Ve)else matrix(fit$Ve,1,1),
    residual_correlation=if(any(is_at))calibration$residual_correlation else NA_real_,
    marker_variance=fit$Vu,residual_scale=fit$Ve,converged=TRUE,
    normal_equations_relative_error=relative_error,cohort_keys=unique(records$cohort))
  attr(values,'training_diagnostics') <- diagnostic
  values
}

redesign_asreml_calibrated <- function(state,frame,inverse) {
  calibration <- redesign_calibrated_residual(state)
  R <- calibration$relative_matrix
  initial <- asreml::asreml(fixed=Pheno~1+Year,random=~vm(Ind,inverse),
    residual=~id(Ind):us(Stage),na.action=asreml::na.method(y='include'),data=frame,start.values=TRUE)
  tab <- initial$vparameters.table
  vcc <- cbind(group=seq_len(nrow(tab)),scale=rep(1,nrow(tab)))
  rownames(vcc) <- tab$Component
  idx <- match(c('Ind:Stage!Stage_PT:PT','Ind:Stage!Stage_AT:PT','Ind:Stage!Stage_AT:AT'),tab$Component)
  if(anyNA(idx))stop('Cannot identify residual parameters for calibrated constraints.')
  vcc[idx,1] <- vcc[idx[1],1]
  vcc[idx,2] <- c(1,R['AT','PT'],R['AT','AT'])
  # Start the three linked parameters at a coherent positive definite matrix.
  scale <- tab$Value[idx[1]]
  initial$R.param[['Ind:Stage']][['Stage']]$initial[] <- c(scale,R['AT','PT']*scale,R['AT','AT']*scale)
  fit <- do.call(asreml::asreml,list(fixed=Pheno~1+Year,random=~vm(Ind,inverse),
    residual=~id(Ind):us(Stage),na.action=asreml::na.method(y='include'),data=frame,
    G.param=initial$G.param,R.param=initial$R.param,vcc=vcc))
  attr(fit,'design_calibration') <- calibration
  fit
}

redesign_fit <- function(state,method,candidates) {
  if(!redesign_uses_repeated_Y9(state)||method=='Pheno')return(redesign_fit_PT_original(state,method,candidates))
  if(method=='GS')return(redesign_native_repeated_gs(state,redesign_training(state,candidates)))
  redesign_fit_repeated_pedigree(state,method,candidates)
}

redesign_predict <- function(pop,method,model) {
  if(is.null(pop)||method=='Pheno')return(pop)
  if(method=='GS'&&is.numeric(model)&&!is.null(names(model))) {
    values<-unname(model[as.character(pop@id)])
    if(any(!is.finite(values)))stop('Missing calibrated genomic predictions.')
    pop@ebv<-matrix(values,ncol=1);return(pop)
  }
  redesign_predict_PT_original(pop,method,model)
}

