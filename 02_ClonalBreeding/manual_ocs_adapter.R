# Historical adapter filename; the active implementation uses native optiSel functions.
source(file.path(getOption('almond.input_dir'),'ocs_optisel.R'),local=.GlobalEnv)
redesign_ocs <- function(state,candidates,model) {
 pool<-redesign_predict(redesign_merge_latest(list(state$parents,candidates)),'GS',model)
 stopifnot(isTRUE(all.equal(state$config$ocs_max_contribution,.10)))
 bv<-as.numeric(ebv(pool));K<-ocs_kinship(pool,state$ocs_frequency)
 use_year<-state$year+if(isTRUE(state$config$first_parental_use>=6L)) 1L else 0L
 age<-use_year-state$birth[as.character(pool@id)]
 allowed<-ocs_compatibility(pool,locuscompt_position_vector,self_compatible_allele)
 allowed<-allowed & outer(age>=state$config$minimum_mother_age,age>=state$config$minimum_pollen_age,'&')
 diag(allowed)<-FALSE
 N<-state$config$crosses*state$config$progeny
 solution<-optisel_offspring_profile(K,bv,allowed,N,state$config$ocs_alpha,state$config$ocs_max_contribution)
 active<-which(solution$metrics$counts>0)
 details<-list(implementation='optiSel_native_offspring',plan=solution$plan,
   family_sizes=solution$family_sizes,offspring=N,max_individual_contribution=.10,
   max_contributing_parents=NA_integer_,group_kinship=solution$metrics$kinship,
   bound=solution$bound,minK=solution$minK,maxK=solution$maxK,
   contributions=solution$metrics$contribution,fallback_to_diversity_endpoint=FALSE,
   endpoints_collapsed=abs(solution$maxK-solution$minK)<1e-10,
   custom_contribution_refinement=FALSE,
   mother_ages=age[solution$plan[,1]],father_ages=age[solution$plan[,2]])
 details$optimization_audit<-solution$audit;state$speed_ocs_solver<-details
 list(parents=pool[active],new=sum(!pool@id[active]%in%state$parents@id),
   bound=solution$bound,minK=solution$minK,maxK=solution$maxK,alpha=state$config$ocs_alpha,
   metrics=solution$metrics,weighted_age=sum(solution$metrics$contribution*age),
   contributions=data.frame(id=as.character(pool@id[active]),contribution=solution$metrics$contribution[active]),
   cross_ids=matrix(as.character(pool@id[solution$plan]),ncol=2),
   family_sizes=solution$family_sizes,solver_details=details)
}
