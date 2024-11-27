# Fill breeding pipeline

# Set initial yield trials with unique individuals

# Sample year effects
P = runif(12)

# Breeding program
for(cohort in 1:12) {
  cat("  FillPipeline stage:",cohort,"of 12\n")
  # Stage 1 Crossing block
  F1 = randCross(Parents, nCrosses = nCrosses, nProgeny = nProgeny)
  if(cohort < 12){
    # Stage 2 Germinate the seedlings in the nursery
    Seedlings = F1
  }
  if(cohort < 11){
    # Stage 3 Juvenility
    HPT1 = Seedlings
  }
  if(cohort < 10){
    # Stage 4 Juvenility
    HPT2 = HPT1
  }
  if(cohort < 9){
    # Stage 5 Record the HPT yields
    HPT3 = HPT2 
  }
  if(cohort < 8){
    # Stage 6 Record the HPT yields
    
    HPT4 = setPheno(HPT3, varE = VarE, reps = repHPT, p = P[cohort], h2 = c(0.2))
  }
  if(cohort < 7){
    # Stage 7 Select 50 superior individuals and plant as advanced clonal trials (ACT)
    ACT1 = selectInd(HPT4, nInd = nClonesACT, use = "pheno")
  }
  if(cohort < 6){
    # Stage 8 Record ACT yields
    ACT2 = ACT1
  }
  if(cohort < 5){
    # Stage 9 Record ACT yields
    ACT3 = setPheno(ACT2, varE = VarE, reps = repACT, p = P[cohort+3L], h2 = c(0.4))
  }
  if(cohort < 4){
    # Stage 10 Select 3 superior individuals and plant as elite clonal trials (ECT)
    ECT1 = selectInd(ACT3, nInd = nClonesECT, use = "pheno")
  }
  if(cohort < 3){
    # Stage 11 Record ECT yields
    ECT2 = ECT1
  }
  if(cohort < 2){
    # Stage 12 Record ECT yields
    ECT3 = setPheno(ECT2, varE = VarE, reps = repECT, p = P[cohort+6L],h2 = c(0.5))
  }
}

