# Rank and retain individuals using phenotypes or predictions.
# The early, PT and AT population replacements are explicit in advance_year.R.
# Public function names are retained for existing launchers and saved workflows.

# Return the highest-ranked n individuals using phenotype or predicted breeding value.
redesign_choose <- function(pop,n,method,model=NULL) {
  pop <- redesign_predict(pop,method,model)
  selectInd(pop,nInd=n,use=if(method=='Pheno') 'pheno' else 'ebv',simParam=SP)
}

