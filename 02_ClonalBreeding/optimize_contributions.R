# Build current Ne-based OCS inputs from the almond candidate pool and validate allocation.
# Contributions are optimized natively; no manual repair or alternative optimizer is used.
source(file.path(getOption('almond.input_dir'),'ocs_optisel.R'),local=.GlobalEnv)
optimize_parent_contributions <- function(state,candidates,model) {
 pool<-predict_breeding_values(merge_latest_populations(list(state$parents,candidates)),'GS',model)
 stopifnot(isTRUE(all.equal(state$config$ocs_max_contribution,.10)))
 bv<-as.numeric(ebv(pool));K<-ocs_kinship(pool,state$ocs_frequency)
 use_year<-state$year+if(isTRUE(state$config$first_parental_use>=6L)) 1L else 0L
 age<-use_year-state$birth[as.character(pool@id)]
 biological_allowed<-ocs_compatibility(pool,locuscompt_position_vector,self_compatible_allele)
 h1<-pullMarkerHaplo(pool,markers=locuscompt_position_vector,haplo=1,simParam=SP)
 h2<-pullMarkerHaplo(pool,markers=locuscompt_position_vector,haplo=2,simParam=SP)
 a<-apply(h1,1,paste,collapse='_');b<-apply(h2,1,paste,collapse='_');sc<-paste(as.numeric(self_compatible_allele),collapse='_')
 keys<-ifelse(a<=b,paste(a,b,sep=':'),paste(b,a,sep=':'))
 keys[a==sc|b==sc]<-paste0('SCindividual_',which(a==sc|b==sc))
 si_groups<-match(keys,unique(keys))
 for(ii in split(seq_along(si_groups),si_groups))stopifnot(!any(biological_allowed[ii,ii,drop=FALSE]))
 allowed<-biological_allowed
 allowed<-allowed & outer(age>=state$config$minimum_mother_age,age>=state$config$minimum_pollen_age,'&')
 diag(allowed)<-FALSE
 N<-state$config$crosses*state$config$progeny
 reference_kinship<-mean(ocs_kinship(state$parents,state$ocs_frequency))
 generation_interval<-mean(use_year-state$birth[as.character(state$parents@id)])
 capture<-getOption('almond.decision_capture_dir')
 if(!is.null(capture))saveRDS(list(year=state$year,scenario=state$config$scenario_id,K=K,bv=bv,allowed=allowed,biological_allowed=biological_allowed,si_groups=si_groups,S_allele1=a,S_allele2=b,genotype_ids=pool@id,age=age,N=N,Ne=state$config$ocs_target_Ne,reference_kinship=reference_kinship,generation_interval=generation_interval,seed=.Random.seed),file.path(capture,'current_decision_input.rds'))
 solution<-optisel_offspring_profile(K,bv,allowed,N,state$config$ocs_target_Ne,state$config$ocs_max_contribution,
  reference_kinship=reference_kinship,generation_interval=generation_interval,si_groups=si_groups,age=age,adult_cap=if(isTRUE(state$config$speed_breeding)).032 else .10)
 active<-which(solution$metrics$counts>0)
 details<-list(implementation='optiSel_adultcap32_v12',plan=solution$plan,
   family_sizes=solution$family_sizes,offspring=N,max_individual_contribution=.10,
   group_kinship=solution$metrics$kinship,
   bound=solution$bound,
   contributions=solution$metrics$contribution,
   custom_contribution_refinement=FALSE,
   mother_ages=age[solution$plan[,1]],father_ages=age[solution$plan[,2]])
 details$target_Ne<-state$config$ocs_target_Ne;details$reference_kinship<-reference_kinship;details$generation_interval_proxy<-generation_interval
 details$optimization_audit<-solution$audit;state$speed_ocs_solver<-details
 list(parents=pool[active],new=sum(!pool@id[active]%in%state$parents@id),
   bound=solution$bound,
   metrics=solution$metrics,weighted_age=sum(solution$metrics$contribution*age),
   contributions=data.frame(id=as.character(pool@id[active]),contribution=solution$metrics$contribution[active]),
   cross_ids=matrix(as.character(pool@id[solution$plan]),ncol=2),
   family_sizes=solution$family_sizes,solver_details=details)
}
