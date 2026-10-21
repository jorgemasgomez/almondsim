# Rank and retain individuals using phenotypes or predictions.
# The early, PT and AT population replacements are explicit in advance_year.R.
# Function names describe each operation; all launchers use this shared implementation.

# Return the highest-ranked n individuals using phenotype or predicted breeding value.
select_best_individuals <- function(pop,n,method,model=NULL) {
  pop <- predict_breeding_values(pop,method,model)
  selectInd(pop,nInd=n,use=if(method=='Pheno') 'pheno' else 'ebv',simParam=SP)
}

