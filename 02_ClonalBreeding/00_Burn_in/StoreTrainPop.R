# Save pedPop records and training population for GS from year 35

# ACT5@fixEff <- as.integer(rep(year,nInd(ACT5)))


ACT3@fixEff <- as.integer(rep(year,nInd(ACT3)))



if(year == startRecords) {
  trainPop = ACT3
  pedPop = rbind(data.frame(Ind   = c(Parents@id),
                            Sire  = c(Parents@father),
                            Dam   = c(Parents@mother),
                            Year  = year,
                            Stage = c(rep("Parents",Parents@nInd)),
                            Pheno = c(Parents@pheno),
                            GV = c(Parents@gv)),
                 data.frame(Ind   = c(ACT3@id),
                            Sire  = c(ACT3@father),
                            Dam   = c(ACT3@mother),
                            Year  = year,
                            Stage = c(rep("ACT",ACT3@nInd)),
                            Pheno = c(ACT3@pheno),
                            GV = c(ACT3@gv)))
}

if (year > startRecords & year < nBurnin+1) {
  trainPop = c(trainPop,ACT3)
  pedPop = rbind(data.frame(Ind   = c(Parents@id),
                            Sire  = c(Parents@father),
                            Dam   = c(Parents@mother),
                            Year  = year,
                            Stage = c(rep("Parents",Parents@nInd)),
                            Pheno = c(Parents@pheno),
                            GV = c(Parents@gv))
                 ,pedPop,
                 data.frame(Ind   = c(ACT3@id),
                            Sire  = c(ACT3@father),
                            Dam   = c(ACT3@mother),
                            Year  = year,
                            Stage = c(rep("ACT3",ACT3@nInd)),
                            Pheno = c(ACT3@pheno),
                            GV = c(ACT3@gv)))
 }



# Separate HPT4 history for GS: keep the three most recent cohorts.
if (year >= startRecords) {
  HPT4@fixEff = as.integer(rep(year, nInd(HPT4)))
  if (year == startRecords) {
    trainPopHPT4 = HPT4
  } else {
    trainPopHPT4 = c(trainPopHPT4, HPT4)
  }
  retainedHPT4Years = tail(sort(unique(trainPopHPT4@fixEff)), 3)
  trainPopHPT4 = trainPopHPT4[which(trainPopHPT4@fixEff %in% retainedHPT4Years)]
}
