# Script name: Genomic tea clonal breeding program
#
# Authors: Initially developed by Nelson Lubanga, Gregor Gorjanc, Jon Bancic;
# exanded/polished for this publication by Jon Bancic, Philip Greenspoon,
# Chris Gaynor, Gregor Gorjanc
#
# Date Created: 2023-12-06

# ---- Clean environment and load packages ----
# rm(list = ls())
# install.packages(pkgs = "AlphaSimR")
library(package = "AlphaSimR")

# ---- Load global parameters ----
source(file = "GlobalParameters.R")
scenarioName = "Burn_in"

# ---- Create list to store results from reps ----
results = list()

for(REP in 1:reps){
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

  save.image(paste0("../burn_in_folder/Burnin_", REP, ".RData"))
  # Save results from current replicate
  # results = append(results, list(output))
}

# Save results
# saveRDS(results, file = paste0(scenarioName,REP,".rds"))


# ---- Analyze results ----
# source(file = "ANALYZERESULTS.R")


