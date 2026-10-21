# Select replacement parents or dispatch to OCS. The caller decides whether selection is used now or stored for next year.
# Public function names are retained for existing launchers and saved workflows.

# Retain ranked existing parents and recruit the prescribed number of incoming parents in non-OCS schemes.
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

# Return a parent allocation; ordinary replacement and OCS share this decision interface.
redesign_allocate <- function(state,candidates,method,model,first,selection_year) {
 validate_ocs_settings(state$config)
 old <- state$parents
 predicted <- redesign_predict(candidates,method,model)
 score <- if(method=='Pheno') predicted@pheno else predicted@ebv
 accuracy <- suppressWarnings(cor(as.numeric(gv(predicted)),as.numeric(score)))
 if(isTRUE(state$config$ocs_enabled)) {
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

