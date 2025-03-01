# Script name: Genomic tea clonal breeding program
#
# Authors: Initially developed by Nelson Lubanga, Gregor Gorjanc, Jon Bancic;
# exanded/polished for this publication by Jon Bancic, Philip Greenspoon,
# Chris Gaynor, Gregor Gorjanc
#
# Date Created: 2023-12-06


# install.packages(pkgs = "AlphaSimR")
library(package = "AlphaSimR")

# pipeline<-FALSE

if (pipeline) {
  scenarioName = "GS_early_selection_extreme"
  # ---- Future phase ----
  # Replace three early stages with genomic prediction WHY???
  #rm(HPT1, HPT2, HPT3)
  nClonesHPT = 250
  for(year in (nBurnin+1):(nBurnin+nFuture)) {
    cat("  Working on future year:",year,"\n")
    source(file = "RunModel_GS.R")    # Run pedigree model
    source(file = "UpdateParents.R")  # Pick parents
    source(file = "AdvanceYear_GS.R") # Advance yield trials by a year
    source(file = "StoreTrainPop.R")  # Store training population
    # Report results
    output$meanG[year] = meanG(Seedlings)
    output$varG[year]  = varG(Seedlings)}
  # Save results from current replicate
  results = append(results, list(output))
  
  
} else {
  # ---- Clean environment and load packages ----
  rm(list = ls())
  source(file = "compatible_crosses.R")
  # ---- Load global parameters ----
  source(file = "GlobalParameters.R")
  scenarioName = "GS"
  
  # ---- Create list to store results from reps ----
  results = list()
  
  for(REP in 1:nReps){
    cat("Working on REP:", REP,"\n")
    
    # ---- Create a data frame to track key parameters ----
    output = data.frame(year     = 1:nCycles,
                        rep      = rep(REP, nCycles),
                        scenario = rep(scenarioName, nCycles),
                        meanG    = numeric(nCycles),
                        varG     = numeric(nCycles),
                        accSel   = numeric(nCycles))
    
    # ---- Create initial parents ----
    source(file = "CreateParents.R")
    
    # ---- Fill breeding pipeline with unique individuals from initial parents ----
    source(file = "FillPipeline.R")
    
    # ---- Simulate year effects ----
    P = runif(nCycles)
    
    # ---- Burn-in phase ----
    for(year in 1:nBurnin) {
      cat("  Working on burnin year:",year,"\n")
      source(file = "UpdateParents.R")  # Pick parents
      source(file = "AdvanceYear.R")    # Advance yield trials by a year
      source(file = "StoreTrainPop.R")  # Store training population
      # Report results
      output$meanG[year] = meanG(Seedlings)
      output$varG[year]  = varG(Seedlings)
    }
    
    # ---- Future phase ----
    # Replace three early stages with genomic prediction WHY???
    #rm(HPT1, HPT2, HPT3)
    
    for(year in (nBurnin+1):(nBurnin+nFuture)) {
      cat("  Working on future year:",year,"\n")
      source(file = "RunModel_GS.R")    # Run pedigree model
      source(file = "UpdateParents.R")  # Pick parents
      source(file = "AdvanceYear_GS.R") # Advance yield trials by a year
      source(file = "StoreTrainPop.R")  # Store training population
      # Report results
      output$meanG[year] = meanG(Seedlings)
      output$varG[year]  = varG(Seedlings)
    }
    
    # Save results from current replicate
    results = append(results, list(output))
  }
  
  # Save results
  saveRDS(results, file = paste0(scenarioName,".rds"))
  
  # ---- Analyze results ----
  source(file = "ANALYZERESULTS.R")
  
  
}





