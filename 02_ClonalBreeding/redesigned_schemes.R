# Comparison of selection methods and parental recycling after 2/4/6/8
# phenotyping years. Existing historical schemes and results are preserved.
library(AlphaSimR)

redesign_config <- function() {
  list(reps=4L, workers=4L, threads=6L, seed=20261005L,
       burnin=40L, future=40L, parents=32L, crosses=20L, progeny=25L,
       act=25L, ect=8L, entrants=8L, training_cohorts=3L,
       early_selection=FALSE, seedlings_keep=500L,
       opr_max_diversity_loss=NA_real_, opr_candidate_limit=50L,
       ocs_alpha=NA_real_, ocs_restarts=4L, ocs_sweeps=12L,
       data_dir=Sys.getenv('ALMONDSIM_DATA_DIR',unset=''),
       founder_source='bank150', legacy_founders=30L,
       founder_bank='C:/Users/franc/Desktop/Proyectos/Almond_sim/Haplotipos_2026-10-05_60000',
       stages=c('F1','Seedlings','Juvenile2','Juvenile3','HPT1','HPT2',
                'HPT3','HPT4','ACT1','ACT2','ACT3','ECT1','ECT2','ECT3'),
       h2=c(HPT1=.2,HPT2=.2,HPT3=.4,HPT4=.4,ACT1=.5,ACT2=.5,
            ACT3=.6,ECT1=.6,ECT2=.6,ECT3=.6,Parents=.6),
       parent_stage=c(`8`='ECT1',`6`='ACT2',`4`='HPT4',`2`='HPT2'))
}

redesign_scenarios <- function(root=NULL) {
  expected <- expand.grid(method=c('Pheno','Pedigree','GS'),phenotyping_years=c(8L,6L,4L,2L),
                          stringsAsFactors=FALSE)
  expected$early_selection <- FALSE
  early <- expand.grid(method=c('Pedigree','GS'),phenotyping_years=c(8L,2L),
                       stringsAsFactors=FALSE)
  early$early_selection <- TRUE
  expected <- rbind(expected,early)
  expected$progeny_per_cross <- ifelse(expected$early_selection,75L,25L)
  expected$seedlings_keep <- 500L
  expected$folder_name <- paste0(expected$method,'_',expected$phenotyping_years,
                                ifelse(expected$early_selection,'_Early1500',''))
  expected$scenario_id <- paste0(expected$method,'_Y',expected$phenotyping_years,
                                ifelse(expected$early_selection,'_Early1500',''))
  expected$new_parents_per_year <- 8L
  expected$renewal_percent <- 25
  extra <- expected[rep(which(expected$folder_name=='GS_2_Early1500'),3),,drop=FALSE]
  extra$renewal_percent <- c(50,75,100)
  extra$new_parents_per_year <- as.integer(floor(32L*extra$renewal_percent/100))
  extra$folder_name <- paste0(extra$folder_name,'_R',extra$renewal_percent)
  extra$scenario_id <- paste0(extra$scenario_id,'_R',extra$renewal_percent)
  expected <- rbind(expected,extra)
  expected$ocs_alpha <- NA_real_
  ocs <- expected[rep(which(expected$folder_name=='GS_2_Early1500'),3),,drop=FALSE]
  ocs$ocs_alpha <- c(.05,.25,.6)
  ocs$new_parents_per_year <- NA_integer_
  ocs$renewal_percent <- NA_real_
  labels <- c('Diversity','Balanced','Gain')
  ocs$folder_name <- paste0('GS_2_Early1500_OCS_',labels)
  ocs$scenario_id <- paste0('GS_Y2_Early1500_OCS_',labels)
  expected <- rbind(expected,ocs)
  expected$opr_max_diversity_loss <- NA_real_
  opr <- expected[rep(which(expected$folder_name=='GS_2_Early1500'),3),,drop=FALSE]
  opr$opr_max_diversity_loss <- c(.05,.10,.20)
  opr$new_parents_per_year <- NA_integer_
  opr$renewal_percent <- NA_real_
  opr$folder_name <- paste0('GS_2_Early1500_OPR_',c('Diversity','Balanced','Gain'))
  opr$scenario_id <- paste0('GS_Y2_Early1500_OPR_',c('Diversity','Balanced','Gain'))
  expected <- rbind(expected,opr)
  expected$folder_name <- sprintf('%02d_%s',seq_len(nrow(expected)),expected$folder_name)
  if(is.null(root)) return(expected)
  rows <- lapply(seq_len(nrow(expected)),function(index) {
    folder <- file.path(root,expected$folder_name[index])
    settings <- new.env(parent=baseenv())
    sys.source(file.path(folder,'GlobalParameters.R'),envir=settings)
    config <- redesign_config()
    stage <- config$parent_stage[[as.character(settings$parent_phenotyping_years)]]
    stopifnot(settings$nParents==config$parents,
              settings$selection_method==expected$method[index],
              settings$parent_phenotyping_years==expected$phenotyping_years[index],
              settings$parent_selection_stage==stage,
              settings$parent_candidate_age==match(stage,config$stages)-1L,
              settings$parent_mean_h2==config$h2[[stage]],
              identical(settings$early_selection,expected$early_selection[index]),
              settings$progeny_per_cross==expected$progeny_per_cross[index],
              settings$seedlings_keep==expected$seedlings_keep[index],
              identical(settings$new_parents_per_year,expected$new_parents_per_year[index]),
              identical(settings$renewal_percent,expected$renewal_percent[index]),
              identical(settings$ocs_alpha,expected$ocs_alpha[index]))
    if(!is.na(expected$opr_max_diversity_loss[index])) stopifnot(identical(settings$opr_max_diversity_loss,expected$opr_max_diversity_loss[index]), settings$opr_candidate_limit==50L)
    row <- expected[index,,drop=FALSE]
    row$folder <- folder
    row
  })
  do.call(rbind,rows)
}

# Nested means: covariance between successive aggregate errors equals the
# smaller error variance. Equal target h2 retains the previous aggregate.
# Initial h2 is calibrated against founder additive variance. Error variance
# remains fixed thereafter; realized h2 can decline with genetic variance.
redesign_error_step <- function(previous,old_variance,new_variance) {
  stopifnot(new_variance <= old_variance + 1e-8)
  ratio <- new_variance/old_variance
  ratio*previous + rnorm(length(previous),sd=sqrt(max(0,new_variance*(1-ratio))))
}

redesign_measure <- function(pop,stage,state) {
  if (is.null(pop)) return(NULL)
  variance <- state$reference_variance*(1/state$config$h2[[stage]]-1)
  ids <- as.character(pop@id)
  residuals <- numeric(length(ids))
  for (i in seq_along(ids)) {
    prior <- state$errors[[ids[i]]]
    if (is.null(prior)) {
      error <- rnorm(1,sd=sqrt(variance))
    } else {
      error <- redesign_error_step(prior$error,prior$variance,variance)
    }
    state$errors[[ids[i]]] <- list(error=error,variance=variance)
    residuals[i] <- error
  }
  # Environment-averaged genetic value. Aggregate noise includes remaining
  # environmental uncertainty; no additional single-year GxE draw is added.
  pop@pheno <- matrix(as.numeric(gv(pop))+residuals,ncol=1)
  pop@fixEff <- rep(as.integer(state$year),nInd(pop))
  pop
}

redesign_merge_latest <- function(pops) {
  pops <- Filter(Negate(is.null),pops)
  combined <- mergePops(pops)
  combined[which(!duplicated(as.character(combined@id),fromLast=TRUE))]
}

redesign_training <- function(state,candidates) {
  # One latest aggregate per individual, never several overlapping means.
  redesign_merge_latest(c(state$history,list(state$pop$HPT4,state$pop$ACT3,
                                           state$parents,candidates)))
}

redesign_fit <- function(state,method,candidates) {
  if (method=='Pheno') return(NULL)
  train <- redesign_training(state,candidates)
  if (method=='GS') return(RRBLUP(train,useReps=FALSE,simParam=SP))
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
  inverse <- asreml::ainverse(ped)
  asreml::asreml.options(trace=FALSE)
  model <- asreml::asreml(fixed=Pheno~1+Year,random=~vm(Ind,inverse),
                         residual=~units,na.action=asreml::na.method(y='include'),data=frame)
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
  setNames(as.numeric(coefficients)[keep],coefficient_ids[keep])
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
  for (i in length(state$config$stages):2) {
    stage <- state$config$stages[i]
    source <- state$config$stages[i-1]
    state$pop[[stage]] <- previous[[source]]
  }
  state$pop$F1 <- NULL
  for (stage in names(state$config$h2)) {
    if (stage!='Parents') state$pop[[stage]] <- redesign_measure(state$pop[[stage]],stage,state)
  }
}

redesign_record <- function(state,scenario,rep,accuracy,entrants,parent_stage) {
  stages <- c('Parents','ECT3','ECT2','ECT1','ACT3','ACT2','ACT1',
              'HPT4','HPT3','HPT2','HPT1','Seedlings','F1')
  row <- list(year=state$year,rep=rep,scenario=scenario,
              meanG=meanG(state$pop$Seedlings),varG=varG(state$pop$Seedlings),
              accSel=accuracy,newParents=entrants,parentSelectionStage=parent_stage,
              accParents=state$parent_accuracy)
  for (stage in stages) {
    pop <- if(stage=='Parents') state$parents else state$pop[[stage]]
    row[[paste0('meanG_',stage)]] <- meanG(pop)
    row[[paste0('varG_',stage)]] <- varG(pop)
    frequencies <- table(c(pop@mother,pop@father))/(2*nInd(pop))
    row[[paste0('HHI',stage)]] <- sum(frequencies^2)
  }
  haplo <- pullMarkerHaplo(state$pop$Seedlings,markers=locuscompt_position_vector,haplo='all')
  row$allelesSI <- length(unique(apply(haplo,1,paste,collapse='_')))
  dosage <- pullSnpGeno(state$pop$Seedlings,chr=6,simParam=SP)
  frequency <- colMeans(dosage)/2
  row$He_chr6 <- mean(2*frequency*(1-frequency))
  as.data.frame(row,stringsAsFactors=FALSE)
}

redesign_year <- function(state,method,phenotyping_years,rep,record=TRUE,renew=TRUE) {
  redesign_advance(state)
  # Existing parents continue being evaluated as they age, including parents
  # selected from HPT2. They do not retain a two-year mean forever.
  age <- state$year-state$birth[as.character(state$parents@id)]
  counts <- pmin(8L,pmax(1L,age-3L))
  parent_parts <- lapply(sort(unique(counts)),function(count) {
    pop <- state$parents[which(counts==count)]
    stage <- names(state$config$h2)[count]
    redesign_measure(pop,stage,state)
  })
  state$parents <- mergePops(parent_parts)
  target <- state$config$parent_stage[[as.character(phenotyping_years)]]
  candidates <- state$pop[[target]]
  model <- if(!is.null(candidates)) redesign_fit(state,method,candidates) else NULL
  accuracy <- NA_real_
  if (!is.null(state$pop$Seedlings) && method!='Pheno' && !is.null(model)) {
    predicted <- redesign_predict(state$pop$Seedlings,method,model)
    accuracy <- cor(as.numeric(gv(predicted)),as.numeric(ebv(predicted)))
  }
  seedlings_before <- if(is.null(state$pop$Seedlings)) 0L else nInd(state$pop$Seedlings)
  if (isTRUE(state$config$early_selection) && seedlings_before>0L) {
    if(method=='Pheno' || is.null(model)) stop('Early selection requires pedigree or genomic predictions.')
    state$pop$Seedlings <- redesign_choose(state$pop$Seedlings,
                                          min(state$config$seedlings_keep,seedlings_before),method,model)
  }
  state$parent_accuracy <- NA_real_
  if (!is.null(candidates)) {
    predicted <- redesign_predict(candidates,method,model)
    score <- if(method=='Pheno') predicted@pheno else predicted@ebv
    state$parent_accuracy <- cor(as.numeric(gv(predicted)),as.numeric(score))
  }
  # Growth/phenotyping precedes selection and crossing in the same calendar year.
  entrants <- 0L
  ocs_solution <- NULL
  opr_solution <- NULL
  old_parent_count <- nInd(state$parents)
  if (renew && !is.na(state$config$ocs_alpha)) {
    if(method!='GS' || is.null(candidates)) stop('OCS requires GS and an eligible phenotyped cohort.')
    ocs_solution <- redesign_ocs(state,candidates,model)
    state$parents <- ocs_solution$parents
    entrants <- ocs_solution$new
    state$ocs_trace[[as.character(state$year)]] <- ocs_solution[c('bound','minK','maxK','alpha','contributions','cross_ids','weighted_age')]
  } else if (renew && !is.null(state$config$opr_max_diversity_loss) && !is.na(state$config$opr_max_diversity_loss)) {
    if(method!='GS' || is.null(candidates)) stop('OPR requires GS and an eligible cohort.')
    opr_solution <- redesign_opr(state,candidates,model)
    state$parents <- opr_solution$parents
    entrants <- opr_solution$new
  } else if (renew && !is.null(candidates)) {
    updated <- redesign_renew(state$parents,candidates,method,model,state$config)
    state$parents <- updated$parents
    entrants <- updated$new
  }
  if (!is.null(state$pop$ACT1)) state$pop$ACT1 <- redesign_choose(state$pop$ACT1,state$config$act,method,model)
  if (!is.null(state$pop$ECT1)) state$pop$ECT1 <- redesign_choose(state$pop$ECT1,state$config$ect,method,model)
  if(is.null(ocs_solution)) {
    state$pop$F1 <- randCrossGamSI(prms=list(SIPos=locuscompt_position_vector,random_failure=FALSE),
                                 pop=state$parents,nCrosses=state$config$crosses,nProgeny=state$config$progeny)
  } else {
    families <- lapply(seq_len(nrow(ocs_solution$cross_ids)),function(i) {
      family <- randCrossGamSI(prms=list(SIPos=locuscompt_position_vector,random_failure=FALSE),
                              pop=state$parents,nCrosses=1L,nProgeny=state$config$progeny,
                              cross_plan=ocs_solution$cross_ids[i,,drop=FALSE])
      if(!is(family,'Pop') || nInd(family)!=state$config$progeny) stop('OCS compatible family did not meet its planned size.')
      family
    })
    state$pop$F1 <- mergePops(families)
    realized <- table(factor(c(state$pop$F1@mother,state$pop$F1@father),levels=state$parents@id))/(2*nInd(state$pop$F1))
    expected <- ocs_solution$contributions$contribution[match(names(realized),ocs_solution$contributions$id)]
    if(anyNA(expected) || max(abs(as.numeric(realized)-expected))>1e-10) stop('Realized OCS contributions differ from the validated mating plan.')
  }
  state$birth[as.character(state$pop$F1@id)] <- state$year
  if(nInd(state$pop$F1)!=state$config$crosses*state$config$progeny) stop('Incomplete compatible progeny cohort.')
  if (!is.null(state$pop$HPT4)) {
    state$history <- c(state$history,list(state$pop$HPT4))
    state$history <- tail(state$history,state$config$training_cohorts)
  }
  if(record) {
    row <- redesign_record(state,paste0(method,'_Y',phenotyping_years,
                          if(isTRUE(state$config$early_selection)) '_Early1500' else ''),rep,accuracy,entrants,target)
    row$targetNewParents <- state$config$entrants
    row$targetRenewalPercent <- 100*state$config$entrants/state$config$parents
    row$renewalRule <- if(!is.null(opr_solution)) 'OPR_free' else if(is.null(ocs_solution)) 'exact' else 'OCS_free'
    row$oprMaxDiversityLoss <- if(is.null(state$config$opr_max_diversity_loss)) NA_real_ else state$config$opr_max_diversity_loss
    for(key in c('beforeHe','afterHe','minHe','lossPercent','expectedBV','tested','rejectedDiversity')) row[[paste0('opr_',key)]] <- if(is.null(opr_solution)) NA_real_ else opr_solution[[key]]
    row$ocsAlpha <- state$config$ocs_alpha
    row$ocsBound <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$bound
    row$ocsGenomicCoancestry <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$metrics$kinship
    row$ocsMateCoancestry <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$metrics$mate_kinship
    row$ocsExpectedBV <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$metrics$gain
    row$ocsEffectiveParents <- if(is.null(ocs_solution)) NA_real_ else 1/sum(ocs_solution$metrics$contribution^2)
    row$weightedParentalAge <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$weighted_age
    row$parentsUsed <- nInd(state$parents)
    if(!is.null(ocs_solution)) {
      row$targetNewParents <- NA_integer_
      row$targetRenewalPercent <- NA_real_
    }
    row$realizedRenewalPercent <- 100*entrants/nInd(state$parents)
    if(!is.null(state$config$scenario_id)) row$scenario <- state$config$scenario_id
    row$earlySelection <- isTRUE(state$config$early_selection)
    row$progenyPerCross <- state$config$progeny
    row$seedlingsBeforeSelection <- seedlings_before
    row$seedlingsAfterSelection <- nInd(state$pop$Seedlings)
    row$parentPhenotypingYears <- phenotyping_years
    row$parentCandidateAge <- match(target,state$config$stages)-1L
    row$meanParentalAge <- mean(state$year-state$birth[as.character(state$parents@id)])
    row
  } else NULL
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
  missing <- state$config$parents-nInd(state$parents)
  if(missing<0L) stop('Founder population exceeds configured parental capacity.')
  if(missing>0L) {
    mature <- redesign_merge_latest(list(state$pop$ECT1,state$pop$ECT2,state$pop$ECT3))
    eligible <- which(!mature@id %in% state$parents@id)
    if(length(eligible)<missing) stop('Insufficient mature candidates to fill the temporary parental pool.')
    added <- redesign_choose(mature[eligible],missing,'Pheno')
    age <- state$year-state$birth[as.character(added@id)]
    if(any(age<11L)) stop('Temporary additional parents are too young.')
    state$parents <- redesign_merge_latest(list(state$parents,added))
  }
  stopifnot(nInd(state$parents)==state$config$parents)
  invisible(state)
}

redesign_setup <- function(root,config) {
  # Reuse validated founder import and S-locus assignment from the original program.
  source(file.path(root,'00_Burn_in','GlobalParameters.R'),local=.GlobalEnv)
  source(file.path(root,'compatible_crosses.R'),local=.GlobalEnv)
  source(file.path(root,'ocs_discrete.R'),local=.GlobalEnv)
  source(file.path(root,'opr_genomic.R'),local=.GlobalEnv)
  if(identical(config$founder_source,'legacy')) {
    bank <- file.path(redesign_data_dir(root,config),'Haplotypes')
    founder_count <- as.integer(config$legacy_founders)
  } else {
    bank <- redesign_founder_bank(config$founder_bank)
    founder_count <- config$parents
  }
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
  state$errors <- new.env(parent=emptyenv())
  state$parents <- Parents
  state$birth <- setNames(rep(-30L,nInd(Parents)),as.character(Parents@id))
  state$pop <- setNames(vector('list',length(config$stages)),config$stages)
  state$history <- list()
  state$year <- -14L
  state$parents <- redesign_measure(state$parents,'Parents',state)
  for(i in seq_len(14)) {
    state$year <- i-14L
    redesign_year(state,'Pheno',8L,0L,record=FALSE,renew=FALSE)
  }
  redesign_fill_parent_pool(state)
  state
}

redesign_worker <- function(job) {
  exclusion_file <- file.path(dirname(job$directory),'excluded_replicates.txt')
  excluded <- if(file.exists(exclusion_file)) as.integer(readLines(exclusion_file,warn=FALSE)) else integer()
  if(job$rep %in% excluded) {
    dir.create(job$directory,recursive=TRUE,showWarnings=FALSE)
    jsonlite::write_json(list(status='excluded',rep=job$rep,phase=job$phase,
      reason='User reduced the active experiment to four replicates before replicate 5 started.',
      updated=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')),
      file.path(job$directory,paste0('progress_',job$phase,'.json')),pretty=TRUE,auto_unbox=TRUE)
    return(list(ok=TRUE,rep=job$rep,result=NULL,excluded=TRUE))
  }
  .libPaths(job$libraries)
  options(almond.input_dir=job$root,almond.nThreads=job$config$threads)
  source(file.path(job$root,'redesigned_schemes.R'),local=.GlobalEnv)
  source(file.path(job$root,'ocs_discrete.R'),local=.GlobalEnv)
  source(file.path(job$root,'opr_genomic.R'),local=.GlobalEnv)
  source(file.path(job$root,'compatible_crosses.R'),local=.GlobalEnv)
  assign('.Random.seed',job$seed,envir=.GlobalEnv)
  dir.create(job$directory,recursive=TRUE,showWarnings=FALSE)
  timing_file <- file.path(job$directory,paste0('timings_',job$phase,'.tsv'))
  options(almond.timing_file=timing_file,almond.timing_rep=job$rep,
          almond.timing_scenario='initialization',almond.timing_year=0L)
  for(name in c('redesign_setup','redesign_advance','redesign_measure','redesign_fit',
               'redesign_predict','redesign_choose','redesign_renew','redesign_ocs',
               'redesign_opr','randCrossGamSI','redesign_record','redesign_year')) {
    assign(name,redesign_timer(get(name,envir=.GlobalEnv),name),envir=.GlobalEnv)
  }
  options(almond.diagnostic_file=file.path(job$directory,paste0('diagnostics_',job$phase,'.log')))
  status_file <- file.path(job$directory,paste0('progress_',job$phase,'.json'))
  progress <- function(status,scenario,year=0L,error=NULL) {
    options(almond.timing_scenario=scenario,almond.timing_year=year)
    jsonlite::write_json(list(status=status,rep=job$rep,phase=job$phase,
      scenario=scenario,year=year,memory=redesign_memory(),updated=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),error=error),
      status_file,pretty=TRUE,auto_unbox=TRUE)
  }
  progress('running','initialization')
  log <- file(file.path(job$directory,paste0('execution_',job$phase,'.log')),'wt')
  sink(log); on.exit({sink();close(log)},add=TRUE)
  tryCatch({
    base_file <- file.path(job$directory,'burnin_redesigned.RData')
    if(job$phase=='parallel') {
      state <- redesign_setup(job$root,job$config)
      for(name in c('redesign_ocs','redesign_opr','randCrossGamSI')) assign(name,redesign_timer(get(name,envir=.GlobalEnv),name),envir=.GlobalEnv)
      burnin <- vector('list',job$config$burnin)
      for(year in seq_len(job$config$burnin)) {
        state$year <- year
        progress('running','Burn_in',year)
        cat('Burn-in year',year,'\n')
        burnin[[year]] <- redesign_year(state,'Pheno',8L,job$rep)
        burnin[[year]]$scenario <- 'Burn_in'
      }
      save(state,burnin,SP,locuscompt_position_vector,self_compatible_allele,file=base_file)
    }
    outputs <- list()
    selected <- job$scenarios[job$scenarios$method %in% job$methods,,drop=FALSE]
    for(index in seq_len(nrow(selected))) {
      load(base_file,envir=.GlobalEnv)
      state <- get('state',.GlobalEnv)
      # Paired schemes start with the same burn-in and future RNG stream.
      assign('.Random.seed',parallel::nextRNGStream(job$seed),envir=.GlobalEnv)
      method <- selected$method[index]; years <- selected$phenotyping_years[index]
      state$config$ocs_alpha <- selected$ocs_alpha[index]
      state$config$opr_max_diversity_loss <- if(is.null(selected$opr_max_diversity_loss)) NA_real_ else selected$opr_max_diversity_loss[index]
      state$config$early_selection <- selected$early_selection[index]
      state$config$progeny <- selected$progeny_per_cross[index]
      state$config$seedlings_keep <- selected$seedlings_keep[index]
      if(is.na(selected$renewal_percent[index])) {
        state$config$entrants <- selected$new_parents_per_year[index]
      } else {
        state$config$entrants <- as.integer(floor(state$config$parents*selected$renewal_percent[index]/100))
      }
      scenario <- selected$scenario_id[index]
      state$config$scenario_id <- scenario
      cat('Scenario',scenario,'\n')
      records <- vector('list',job$config$future)
      for(k in seq_len(job$config$future)) {
        state$year <- job$config$burnin+k
        progress('running',scenario,state$year)
        cat('Future year',state$year,'\n')
        records[[k]] <- redesign_year(state,method,years,job$rep)
      }
      outputs[[scenario]] <- do.call(rbind,records)
      saveRDS(outputs[[scenario]],file.path(job$directory,paste0('results_',scenario,'.rds')))
      progress('scenario_complete',scenario,state$year)
      if(!is.na(state$config$ocs_alpha)) saveRDS(state$ocs_trace,file.path(job$directory,paste0('ocs_plan_',scenario,'.rds')))
    }
    result <- do.call(rbind,outputs)
    saveRDS(result,file.path(job$directory,paste0('results_',job$phase,'.rds')))
    progress('complete','all')
    list(ok=TRUE,rep=job$rep,result=result)
  },error=function(e) {
    progress('failed',getOption('almond.timing_scenario'),getOption('almond.timing_year'),conditionMessage(e))
    writeLines(conditionMessage(e),file.path(job$directory,paste0('error_',job$phase,'.txt')))
    list(ok=FALSE,rep=job$rep,error=conditionMessage(e))
  })
}

run_redesigned_schemes <- function(root,config=redesign_config(),scenarios=redesign_scenarios()) {
  experiment_start <- proc.time()[['elapsed']]
  root <- normalizePath(root,winslash='/',mustWork=TRUE)
  if(any(scenarios$method=='Pedigree')) redesign_check_asreml()
  if(!identical(config$founder_source,'legacy')) redesign_founder_bank(config$founder_bank)
  dir.create(file.path(redesign_data_dir(root,config),'redesigned_runs'),recursive=TRUE,showWarnings=FALSE)
  run <- tempfile(paste0('run_',format(Sys.time(),'%Y%m%d_%H%M%S'),'_'),
                  tmpdir=file.path(redesign_data_dir(root,config),'redesigned_runs'))
  dir.create(run,recursive=TRUE)
  writeLines(run,file.path(redesign_data_dir(root,config),'latest_run.txt'))
  RNGkind("L'Ecuyer-CMRG");set.seed(config$seed)
  stream <- .Random.seed
  jobs <- vector('list',config$reps)
  for(rep in seq_len(config$reps)) {
    jobs[[rep]] <- list(rep=rep,root=root,config=config,scenarios=scenarios,
                       methods=c('Pheno','GS'),phase='parallel',seed=stream,
                       libraries=.libPaths(),directory=file.path(run,sprintf('rep_%03d',rep)))
    stream <- parallel::nextRNGStream(stream)
  }
  saveRDS(list(config=config,scenarios=scenarios),file.path(run,'configuration.rds'))
  write.csv(scenarios,file.path(run,'scenarios.csv'),row.names=FALSE)
  cluster <- parallel::makeCluster(min(config$workers,config$reps),outfile=file.path(run,'cluster.log'))
  on.exit(parallel::stopCluster(cluster),add=TRUE)
  results <- parallel::parLapplyLB(cluster,jobs,function(job) {
    .libPaths(job$libraries)
    source(file.path(job$root,'redesigned_schemes.R'))
    redesign_worker(job)
  })
  if(any(!vapply(results,`[[`,logical(1),'ok'))) stop('A Pheno/GS job failed; inspect ',run)
  # All pedigree schedules run sequentially in one licensed R worker.
  pedigree <- vector('list',config$reps)
  for(rep in seq_len(config$reps)) {
    job <- jobs[[rep]];job$methods <- 'Pedigree';job$phase <- 'pedigree'
    pedigree[[rep]] <- parallel::parLapply(cluster[1],list(job),function(job) {
      .libPaths(job$libraries)
    source(file.path(job$root,'redesigned_schemes.R'))
      redesign_worker(job)
    })[[1]]
    if(!pedigree[[rep]]$ok) stop('Pedigree failed: ',pedigree[[rep]]$error,'; ',run)
  }
  combined <- do.call(rbind,lapply(c(results,pedigree),`[[`,'result'))
  saveRDS(combined,file.path(run,'results_all_schemes.rds'))
  write.table(combined,file.path(run,'results_all_schemes.txt'),sep='\t',row.names=FALSE,quote=FALSE)
  saveRDS(list(complete=TRUE,replicates=length(unique(combined$rep)),scenarios=scenarios),file.path(run,'status.rds'))
  saveRDS(list(elapsed_seconds=proc.time()[['elapsed']]-experiment_start,completed=Sys.time()),file.path(run,'total_runtime.rds'))
  message('Complete redesigned experiment: ',run)
  invisible(list(run_dir=run,data=combined))
}



redesign_data_dir <- function(root,config=list()) {
  configured <- config$data_dir
  if(is.null(configured) || !nzchar(configured)) configured <- Sys.getenv('ALMONDSIM_DATA_DIR',unset='')
  if(!nzchar(configured)) configured <- file.path(dirname(dirname(normalizePath(root,winslash='/',mustWork=TRUE))),'almondsim_data')
  normalizePath(configured,winslash='/',mustWork=FALSE)
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
