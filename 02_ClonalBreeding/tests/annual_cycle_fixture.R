# Small deterministic populations for testing annual orchestration only.
# Prediction, phenotyping, crossing and parent allocation are controlled test
# doubles. Their numerical algorithms have their own checks; this fixture tests
# call order, population replacement, parental timing, output fields and RNG use.
if (!methods::isClass('AnnualCyclePopulation')) {
  methods::setClass('AnnualCyclePopulation', slots=c(
    id='character', pheno='matrix', ebv='matrix', truth='numeric',
    mother='character', father='character', misc='list'))
  methods::setMethod('[',signature(x='AnnualCyclePopulation',i='ANY',j='missing',drop='missing'),
    function(x,i,j,...,drop)methods::new('AnnualCyclePopulation',
      id=x@id[i],pheno=x@pheno[i,,drop=FALSE],ebv=x@ebv[i,,drop=FALSE],
      truth=x@truth[i],mother=x@mother[i],father=x@father[i],misc=x@misc))
}

annual_fixture_population <- function(ids, evaluated=TRUE, truth=NULL) {
  if(is.null(truth))truth <- rep(c(4,1,3,2),length.out=length(ids))
  methods::new('AnnualCyclePopulation',id=ids,truth=as.numeric(truth),
    pheno=matrix(if(evaluated)truth else rep(NA_real_,length(ids))),
    ebv=matrix(truth),mother=rep('',length(ids)),father=rep('',length(ids)),
    misc=list(training_record_stage=rep('PT',length(ids))))
}

annual_fixture_hash <- function(value) {
  path <- tempfile();on.exit(unlink(path))
  saveRDS(value,path,version=2,compress=FALSE)
  unname(tools::md5sum(path))
}

annual_fixture_engine <- function(root, year_function=NULL) {
  engine <- new.env(parent=.GlobalEnv)
  sys.source(file.path(root,'redesigned_schemes.R'),envir=engine)
  if(!is.null(year_function)) {
    environment(year_function)<-engine;engine$advance_breeding_year<-year_function
  }
  engine$trace <- character()
  log <- function(text)engine$trace<-c(engine$trace,text)
  engine$nInd <- function(pop)length(pop@id)
  engine$gv <- function(pop)pop@truth
  engine$ebv <- function(pop)pop@ebv
  engine$meanG <- function(pop)mean(pop@truth)
  engine$varG <- function(pop)stats::var(pop@truth)
  engine$pullMarkerHaplo <- function(pop,...)matrix(seq_along(pop@id),ncol=1)
  engine$locuscompt_position_vector <- 1L
  engine$SP <- NULL

  engine$predict_breeding_values <- function(pop,method,model) {
    log(paste('predict',method,pop@id[1],if(is.null(model))'none' else model$label))
    if(method!='Pheno')pop@ebv<-matrix(pop@truth+if(model$label=='AT3')0.2 else 0.1)
    pop
  }
  engine$selectInd <- function(pop,nInd,use,simParam) {
    log(paste('select',pop@id[1],nInd,use))
    score<-if(use=='pheno')pop@pheno else pop@ebv
    pop[order(as.numeric(score),decreasing=TRUE)[seq_len(nInd)]]
  }
  engine$fit_decision_model <- function(state,method,anchor) {
    log(paste('fit',method,anchor));if(method=='Pheno')NULL else list(label=anchor)
  }
  engine$validate_decision_training <- function(pop,model,anchor) {
    log(paste('check_training',anchor,pop@id[1]))
    stopifnot(identical(model$label,anchor));invisible(TRUE)
  }
  engine$update_phenotypes <- function(pop,stage,state) {
    log(paste('measure',stage));if(is.null(pop))return(NULL)
    pop@pheno<-matrix(pop@truth+state$year/100+rnorm(length(pop@id),sd=.01));pop
  }
  engine$update_parent_phenotypes <- function(state) {
    state$parents<-engine$update_phenotypes(state$parents,'Parents',state);invisible(state)
  }
  engine$archive_evaluated_cohorts <- function(state) {
    log('archive')
    for(stage in c('PT2','AT3'))state$stage_history[[stage]][[as.character(state$year)]]<-state$pop[[stage]]
    invisible(state)
  }
  engine$build_training_population <- function(state,candidates,anchor=NULL) {
    log('training_metadata');annual_fixture_population(paste0('training',1:6))
  }
  make_allocation <- function(state,parents,first,selection_year) {
    contribution<-c(.5,.5,rep(0,length(parents@id)-2L))
    ocs<-if(!state$config$ocs_enabled)NULL else list(
      parents=parents,new=3L,bound=.2,metrics=list(kinship=.1,mate_kinship=.05,
        gain=mean(parents@truth),contribution=contribution),
      contributions=data.frame(id=parents@id,contribution=contribution),
      cross_ids=matrix(parents@id[1:2],ncol=2),family_sizes=4L,
      solver_details=list(optimization_audit=list(rounded_kinship_excess_relative=0)))
    list(parents=parents,new=3L,ocs=ocs,selected_year=selection_year,
         use_year=selection_year+if(first>=6L)1L else 0L,accuracy=.7)
  }
  engine$fixture_allocation <- make_allocation
  engine$select_parent_allocation <- function(state,candidates,method,model,first,selection_year) {
    log(paste('parents',first,if(state$config$ocs_enabled)'OCS' else 'replacement'))
    make_allocation(state,candidates[1:3],first,selection_year)
  }
  engine$randCrossGamSI <- function(settings,parents,nCrosses,nProgeny,cross_plan=NULL,family_sizes=NULL) {
    state<-getOption('almond.cross_state');log('cross')
    n<-if(is.null(family_sizes))nCrosses*nProgeny else sum(family_sizes)
    pop<-annual_fixture_population(paste0('offspring_',state$year,'_',seq_len(n)),FALSE,runif(n))
    pop@mother<-rep(parents@id[1],n);pop@father<-rep(parents@id[2],n)
    state$speed_cross_trace<-list(distinct_pair_count=1L,
      mother_ages=rep(state$year-state$birth[parents@id[1]],n),
      father_ages=rep(state$year-state$birth[parents@id[2]],n))
    pop
  }
  engine
}

annual_fixture_case <- function(engine,scenario,cached=TRUE,record=TRUE,renew=TRUE) {
  state<-new.env(parent=emptyenv());config<-engine$breeding_config()
  fields<-setdiff(intersect(names(config),names(scenario)),c('training_stage','training_cohort_stage'))
  for(field in fields)config[[field]]<-scenario[[field]]
  config$scenario_id<-scenario$scenario_id
  config$first_parental_use<-scenario$first_parental_use
  config$parents<-3L;config$crosses<-2L;config$progeny<-2L
  config$seedlings_keep<-3L;config$act<-2L;config$ect<-1L
  config$minimum_pollen_age<-if(scenario$speed_breeding)2L else 4L
  state$config<-config;state$year<-41L
  state$parents<-annual_fixture_population(paste0('parent',1:3))
  state$pop<-setNames(lapply(seq_along(config$stages),function(i)
    annual_fixture_population(paste0(config$stages[i],'_',1:4),i>=5L)),config$stages)
  state$birth<-setNames(rep(20L,3),state$parents@id)
  for(i in seq_along(config$stages))state$birth[state$pop[[i]]@id]<-40L-(i-1L)
  state$environment_p<-c(`41`=.25)
  state$stage_history<-list(PT2=as.list(1:3),AT3=as.list(1:3))
  state$ocs_trace<-list();state$pending<-NULL
  if(renew&&scenario$first_parental_use>=6L) {
    parents<-annual_fixture_population(paste0('pending',1:3));state$birth[parents@id]<-30L
    state$pending<-engine$fixture_allocation(state,parents,scenario$first_parental_use,40L)
  }
  state$model_method<-if(cached)scenario$method else NULL
  state$prediction_pt_model<-if(cached)list(label='previous_PT2') else NULL
  engine$trace<-character();set.seed(84321L)
  row<-engine$advance_breeding_year(state,scenario$method,scenario$first_parental_use,1L,record,renew)
  # This fixture controls phenotyping with a test double; compare orchestration
  # with the original reference, normalizing only the explicit model labels.
  # Real stage phenotyping is checked in check_stage_phenotyping.R.
  if(!is.null(row))row$phenotypingModel<-NULL
  state$config$phenotyping_model<-NULL
  state$config$design_version<-'PT_AT_ET_optiSel_adultcap32_v12'
  if(!is.null(row))row$designVersion<-'PT_AT_ET_optiSel_adultcap32_v12'
  snapshot<-as.list(state);snapshot<-snapshot[sort(names(snapshot))]
  c(row_md5=annual_fixture_hash(row),state_md5=annual_fixture_hash(snapshot),
    calls_md5=annual_fixture_hash(engine$trace),rng_md5=annual_fixture_hash(.Random.seed))
}

annual_fixture_results <- function(engine,root) {
  catalogue<-engine$all_breeding_scenarios(root);results<-list()
  for(i in seq_len(nrow(catalogue)))for(cached in c(TRUE,FALSE)) {
    key<-paste(catalogue$scenario_id[i],if(cached)'cached' else 'initial',sep=':')
    results[[key]]<-annual_fixture_case(engine,as.list(catalogue[i,]),cached)
  }
  baseline<-as.list(catalogue[1,])
  for(record in c(TRUE,FALSE))for(renew in c(TRUE,FALSE)) {
    key<-paste('initialization',record,renew,sep=':')
    results[[key]]<-annual_fixture_case(engine,baseline,FALSE,record,renew)
  }
  data.frame(case=names(results),do.call(rbind,results),row.names=NULL,check.names=FALSE)
}
