# Hypothetical pollen production at age 2; seed production remains age >=4.
speed_allowed <- function(compatibility,age) {
 stopifnot(length(age)==nrow(compatibility),!anyNA(age))
 allowed<-compatibility & outer(age>=4L,age>=2L,'&')
 diag(allowed)<-FALSE
 allowed
}
speed_measure_parents <- function(state) {
 age<-state$year-state$birth[as.character(state$parents@id)]
 if(anyNA(age)||any(age<2L)) stop('Unknown or ineligible parental age.')
 parts<-list()
 juvenile<-which(age<4L)
 if(length(juvenile)) {
   pop<-state$parents[juvenile]
   pop@pheno[,]<-NA_real_
   parts[[1L]]<-pop
 }
 counts<-pmin(8L,age-3L)
 for(count in sort(unique(counts[age>=4L]))) {
   pop<-state$parents[which(age>=4L & counts==count)]
   parts[[length(parts)+1L]]<-redesign_measure(pop,names(state$config$h2)[count],state)
 }
 state$parents<-mergePops(parts)
}
speed_install <- function(adapter) {
 source(adapter,local=.GlobalEnv)
 source(file.path(dirname(adapter),'optisel_roles.R'),local=.GlobalEnv)
 original_training<-redesign_training
 assign('redesign_training',function(state,candidates) {
   train<-original_training(state,candidates)
   keep<-which(is.finite(as.numeric(train@pheno)))
   if(!length(keep)) stop('No phenotyped individuals for genomic training.')
   train[keep]
 },.GlobalEnv)
 original_cross<-randCrossGamSI
 assign('randCrossGamSI',function(prms,pop,nCrosses=NA,nProgeny,cross_plan=NA) {
   state<-getOption('almond.speed_state')
   age<-state$year-state$birth[as.character(pop@id)]
   allowed<-speed_allowed(ocs_compatibility(pop,locuscompt_position_vector,self_compatible_allele),age)
   if(is.na(cross_plan[1])) {
     edges<-which(allowed,arr.ind=TRUE)
     if(!nrow(edges)) stop('No compatible mature mother / eligible pollen donor pair.')
     plan<-edges[sample.int(nrow(edges),nCrosses,replace=TRUE),,drop=FALSE]
   } else {
     plan<-cross_plan
     if(is.character(plan)) plan<-matrix(match(plan,as.character(pop@id)),ncol=2)
   }
   if(anyNA(plan)||!all(allowed[plan])) stop('Speed breeding plan violates age or S compatibility.')
   families<-lapply(seq_len(nrow(plan)),function(i)
     original_cross(prms,pop,nCrosses=1L,nProgeny=nProgeny,cross_plan=plan[i,,drop=FALSE]))
   out<-mergePops(families)
   mother_age<-state$year-state$birth[as.character(out@mother)]
   father_age<-state$year-state$birth[as.character(out@father)]
   stopifnot(nInd(out)==nrow(plan)*nProgeny,all(mother_age>=4),all(father_age>=2))
   previous<-state$speed_cross_trace
   mother_age<-c(previous$mother_ages,mother_age)
   father_age<-c(previous$father_ages,father_age)
   state$speed_cross_trace<-list(mother_ages=mother_age,father_ages=father_age,
      mother_age=mean(mother_age),father_age=mean(father_age),
      pollen_from_age2=mean(father_age==2),pollen_from_under4=mean(father_age<4),
      cross_ids=rbind(previous$cross_ids,matrix(as.character(pop@id[plan]),ncol=2)))
   out
 },.GlobalEnv)
 # Replacement quota remains exactly 8/32 in the arm without OCS.
 original_renew<-redesign_renew
 assign('redesign_renew',function(parents,candidates,method,model,config) {
   updated<-original_renew(parents,candidates,method,model,config)
   state<-getOption('almond.speed_state')
   age<-state$year-state$birth[as.character(updated$parents@id)]
   if(!any(speed_allowed(ocs_compatibility(updated$parents,locuscompt_position_vector,self_compatible_allele),age)))
     stop('Selected parents have no compatible mature mother; no invalid cross made.')
   updated
 },.GlobalEnv)
 assign('redesign_ocs',function(state,candidates,model) {
   pool<-redesign_predict(redesign_merge_latest(list(state$parents,candidates)),'GS',model)
   bv<-as.numeric(ebv(pool));K<-ocs_kinship(pool,state$ocs_frequency)
   age<-state$year-state$birth[as.character(pool@id)]
   allowed<-ocs_compatibility(pool,locuscompt_position_vector,self_compatible_allele)
   if(isTRUE(state$config$speed_breeding)) allowed<-speed_allowed(allowed,age)
   scaled<-if(sd(bv)>0) (bv-mean(bv))/sd(bv) else rep(0,length(bv))
   maternal<-if(isTRUE(state$config$speed_breeding)) age>=4L else rep(TRUE,length(age))
   paternal<-if(isTRUE(state$config$speed_breeding)) age>=2L else rep(TRUE,length(age))
   solution<-optisel_role_solution(K,scaled,allowed,maternal,paternal,state$config$ocs_alpha,
     state$config$crosses,state$config$parents)
   bound<-solution$bound;minK<-solution$minK;maxK<-solution$maxK
   state$speed_ocs_solver<-solution
   active<-which(solution$contributions>0)
   metrics<-ocs_metrics(solution$plan,K,bv)
   list(parents=pool[active],new=sum(!pool@id[active]%in%state$parents@id),
     bound=bound,minK=minK,maxK=maxK,alpha=state$config$ocs_alpha,metrics=metrics,
     weighted_age=sum(solution$contributions*age),
     contributions=data.frame(id=as.character(pool@id[active]),contribution=solution$contributions[active]),
     cross_ids=matrix(as.character(pool@id[solution$plan]),ncol=2),solver_details=solution)
 },.GlobalEnv)
 # Preserve the existing pipeline and metrics; replace only parental evaluation.
 year_function<-redesign_year
 code<-deparse(body(year_function),width.cutoff=500L)
 start<-grep('^    age <- state\\$year',code)
 end<-grep('^    state\\$parents <- mergePops\\(parent_parts\\)',code)
 stopifnot(length(start)==1L,length(end)==1L,end>start)
 code<-c(code[seq_len(start-1L)],'    speed_measure_parents(state)',code[(end+1L):length(code)])
 body(year_function)<-parse(text=paste(code,collapse='\n'))[[1]]
 assign('redesign_year',function(state,method,phenotyping_years,rep,record=TRUE,renew=TRUE) {
   options(almond.speed_state=state)
   on.exit(options(almond.speed_state=NULL))
   state$speed_cross_trace<-NULL
   row<-year_function(state,method,phenotyping_years,rep,record,renew)
   if(record) {
     row$parentPhenotypingYears<-0L
     row$parentCandidateAge<-2L
     row$minimumPollenAge<-2L;row$minimumSeedParentAge<-4L
     row$meanMotherAge<-state$speed_cross_trace$mother_age
     row$meanFatherAge<-state$speed_cross_trace$father_age
     row$weightedParentalAge<-(row$meanMotherAge+row$meanFatherAge)/2
     row$pollenFromAge2<-state$speed_cross_trace$pollen_from_age2
     row$pollenFromUnder4<-state$speed_cross_trace$pollen_from_under4
     row$ocsImplementation<-if(is.na(state$config$ocs_alpha)) NA_character_ else 'optiSel + lpSolve compatibility recovery'
   }
   row
 },.GlobalEnv)
}

optisel_install <- function(adapter) {
 names<-c('redesign_year','redesign_training','redesign_renew','randCrossGamSI')
 originals<-mget(names,envir=.GlobalEnv)
 speed_install(adapter)
 for(name in names) assign(name,originals[[name]],envir=.GlobalEnv)
}
