# Save pedPop records and training population for GS from year 35

ACT3@fixEff <- as.integer(rep(year,nInd(ACT3)))
# ECT6@fixEff <- as.integer(rep(year,nInd(ECT6)))

if(year == startRecords) {
  # trainPop = ECT6
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
  # trainPop = c(trainPop,ECT6)
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

if (year > nBurnin) {
  # Update training population (set to keep 6 years worth of records)
  # remove   = table(trainPop@fixEff)[1]            # remove oldest records
  # trainPop = c(trainPop[-c(1:remove)], ECT6)      # add new records

  # Update pedPop with new records from ACT stage
  pedPop = rbind(data.frame(Ind   = c(Parents@id),
                            Sire  = c(Parents@father),
                            Dam   = c(Parents@mother),
                            Year  = year,
                            Stage = c(rep("Parents",Parents@nInd)),
                            Pheno = c(Parents@pheno),
                            GV = c(Parents@gv)),
                 pedPop,
                 data.frame(Ind   = c(ACT3@id),
                            Sire  = c(ACT3@father),
                            Dam   = c(ACT3@mother),
                            Year  = year,
                            Stage = c(rep("ACT",ACT3@nInd)),
                            Pheno = c(ACT3@pheno),
                            GV = c(ACT3@gv)))
}
