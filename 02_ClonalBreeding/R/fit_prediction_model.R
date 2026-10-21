# Fit phenotypic, pedigree and genomic decision models. Explicit PT-only implementations support calibrated repeated PT/AT models in year-9 schemes.
# Public function names are retained for existing launchers and saved workflows.

# Fit native AlphaSimR RRBLUP to the PT-only training set and attach cohort diagnostics.
redesign_native_rrblup <- function(train,state) {
  if(any(train@misc$training_record_stage!='PT')) stop('GS training must contain only original PT records.')
  model <- AlphaSimR::RRBLUP(train,simParam=SP)
  attr(model,'training_diagnostics') <- list(solver='AlphaSimR::RRBLUP',
    version=as.character(packageVersion('AlphaSimR')),record_counts=table(train@misc$training_record_stage),
    cohort_keys=unique(train@misc$training_record_cohort),fixed_year=TRUE)
  model
}

# Reject a decision model with the wrong anchor or missing candidate own PT records.
redesign_assert_decision_training <- function(pop,model,anchor) {
  if(is.null(pop)) return(invisible(TRUE))
  meta <- attr(model,'decision_training')
  if(is.null(meta) || !identical(meta$anchor,anchor) ||
     !all(as.character(pop@id) %in% meta$ids))
    stop('Decision model excludes own PT records or uses the wrong cohort anchor: ',anchor)
  invisible(TRUE)
}

# Select the PT-only decision anchor without changing model timing.
redesign_fit_for_stage_PT_original <- function(state,method,anchor) {
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

# Fit the PT-only model. Year-9 dispatch may call the calibrated repeated-record implementation instead.
redesign_fit_PT_original <- function(state,method,candidates) {
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

# Provide the explicit validation-only pedigree solver; production pedigree uses ASReml.
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

# Build the pedigree repeated-record fit with the same cohort and residual design as genomic prediction.
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

# Fit and audit the model for one selection anchor, dispatching PT-only or year-9 repeated records.
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

# Return the prespecified PT/AT relative residual covariance from the simulation error design.
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

# Whiten repeated responses and designs, fit solveRRBLUP, and validate its normal equations.
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

# Constrain ASReml's PT/AT residual parameters to a shared scale and fixed relative covariance.
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

# Dispatch to PT-only, repeated genomic, or repeated pedigree prediction without changing candidate timing.
redesign_fit <- function(state,method,candidates) {
  if(!redesign_uses_repeated_Y9(state)||method=='Pheno')return(redesign_fit_PT_original(state,method,candidates))
  if(method=='GS')return(redesign_native_repeated_gs(state,redesign_training(state,candidates)))
  redesign_fit_repeated_pedigree(state,method,candidates)
}

