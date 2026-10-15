source(file.path(getOption('almond.input_dir'),'ocs_optisel_roles.R'),local=.GlobalEnv)
source(file.path(getOption('almond.input_dir'),'ocs_optisel.R'),local=.GlobalEnv)
# optiSel contributions and mating allocation with a hard combined 10% tree cap.
redesign_ocs <- function(state,candidates,model) {
 pool<-redesign_predict(redesign_merge_latest(list(state$parents,candidates)),'GS',model)
 stopifnot(isTRUE(all.equal(state$config$ocs_max_contribution,.10)))
 bv<-as.numeric(ebv(pool));K<-ocs_kinship(pool,state$ocs_frequency)
 use_year<-state$year+if(isTRUE(state$config$first_parental_use>=6L)) 1L else 0L
 age<-use_year-state$birth[as.character(pool@id)]
 allowed<-ocs_compatibility(pool,locuscompt_position_vector,self_compatible_allele)
 allowed<-allowed & outer(age>=state$config$minimum_mother_age,age>=state$config$minimum_pollen_age,'&')
 diag(allowed)<-FALSE
 solution<-optisel_matings_profile(K,bv,allowed,state$config$crosses,state$config$parents,
   state$config$ocs_alpha,restarts=state$config$ocs_restarts,sweeps=state$config$ocs_sweeps,
   seed=state$config$seed+as.integer(state$year))
 active<-which(solution$metrics$counts>0)
 stopifnot(all(allowed[solution$plan]),solution$metrics$kinship<=solution$bound+1e-8)
 details<-list(implementation='optiSel_contributions_and_matings',plan=solution$plan,
   max_individual_contribution=.10,max_participations=ocs_cap_slots(state$config$crosses),
   group_kinship=solution$metrics$kinship,bound=solution$bound,minK=solution$minK,maxK=solution$maxK,
   contributions=solution$metrics$contribution,fallback_to_diversity_endpoint=FALSE,
   endpoints_collapsed=abs(solution$maxK-solution$minK)<1e-10,
   contribution_recovery=FALSE,diversity_endpoint_recovered=FALSE,
   recovery_stage='see_optimization_audit',participations_reassigned=0L,
   attempts=state$config$ocs_restarts,mother_ages=age[solution$plan[,1]],
   father_ages=age[solution$plan[,2]])
 details$optimization_audit<-solution$audit; state$speed_ocs_solver<-details
 list(parents=pool[active],new=sum(!pool@id[active]%in%state$parents@id),
   bound=solution$bound,minK=solution$minK,maxK=solution$maxK,alpha=state$config$ocs_alpha,
   metrics=solution$metrics,weighted_age=sum(solution$metrics$contribution*age),
   contributions=data.frame(id=as.character(pool@id[active]),contribution=solution$metrics$contribution[active]),
   cross_ids=matrix(as.character(pool@id[solution$plan]),ncol=2),solver_details=details)
}
