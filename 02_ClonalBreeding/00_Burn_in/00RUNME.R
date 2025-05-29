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
                      accSel   = numeric(nCycles),
                      HHIact1 = numeric(nCycles),
                      HHIect1 = numeric(nCycles)
                      )

  # ---- Create initial parents ----
  source(file = "CreateParents.R")

  # ---- Fill breeding pipeline with unique individuals from initial parents ----
  source(file = "FillPipeline.R")

  # ---- Simulate year effects ----
  P = runif(nCycles)
  
  # Lista donde guardar HHI por año
  parents_hhi_year_list <- list()
  # lista donde se irán guardando los padres de cada año
  
  # ---- Burn-in phase ----
  for(year in 1:nBurnin) {
    cat("  Working on burnin year:",year,"\n")
    source(file = "UpdateParents.R")  # Pick parents
    source(file = "AdvanceYear.R")    # Advance yield trials by a year
    source(file = "StoreTrainPop.R")  # Store training population
    # Report results
    output$meanG[year] = meanG(Seedlings)
    output$varG[year]  = varG(Seedlings)
    
    
    # Obtener IDs combinados (madre + padre)
    act1_ids <- c(as.character(ACT1@mother), as.character(ACT1@father))
    ect1_ids <- c(as.character(ECT1@mother), as.character(ECT1@father))
    
    
    # Calcular HHI
    act1_freq <- table(act1_ids) / length(act1_ids)
    print(act1_freq)
    act1_hhi_raw <- sum(act1_freq^2)
    act1_hhi_norm <- act1_hhi_raw
    # act1_hhi_norm <- act1_hhi_raw / (1 / (nParents*2))
    
    ect1_freq <- table(ect1_ids) / length(ect1_ids)
    ect1_hhi_raw <- sum(ect1_freq^2)
    ect1_hhi_norm <- ect1_hhi_raw
    # ect1_hhi_norm <- ect1_hhi_raw / (1 / (nParents*2))
    
    
    output$HHIact1[year]=act1_hhi_norm
    output$HHIect1[year]=ect1_hhi_norm
    
  }

  save.image(paste0("../burn_in_folder/Burnin_", REP, ".RData"))
  # Save results from current replicate
  # results = append(results, list(output))
}



# Save results
# saveRDS(results, file = paste0(scenarioName,REP,".rds"))


# ---- Analyze results ----
# source(file = "ANALYZERESULTS.R")


