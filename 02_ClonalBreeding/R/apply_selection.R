# Rank individuals and apply the fixed PT and AT trial cuts. Early selection uses the same ranking helper, called from advance_year.R.
# Public function names are retained for existing launchers and saved workflows.

# Return the highest-ranked n individuals using phenotype or predicted breeding value.
redesign_choose <- function(pop,n,method,model=NULL) {
  pop <- redesign_predict(pop,method,model)
  selectInd(pop,nInd=n,use=if(method=='Pheno') 'pheno' else 'ebv',simParam=SP)
}

# Apply PT2 -> 25 and AT3 -> 8 only after training archives and parental decisions are updated.
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

