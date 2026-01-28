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
  
  stages <- c("Parents","ECT3", "ECT2", "ECT1", "ACT3", "ACT2", "ACT1", 
              "HPT4", "HPT3", "HPT2", "HPT1", "Seedlings", "F1")

  # ---- Create a data frame to track key parameters ----
  # Create all meanG and varG column names for each stage
  meanG_cols <- paste0("meanG_", stages)
  varG_cols  <- paste0("varG_", stages)
  
  output = data.frame(
    year     = 1:nCycles,
    rep      = rep(REP, nCycles),
    scenario = rep(scenarioName, nCycles),
    meanG    = numeric(nCycles),
    varG     = numeric(nCycles),
    accSel   = numeric(nCycles),
    HHIACT1  = numeric(nCycles),
    HHIECT1  = numeric(nCycles),
    HHIACT2  = numeric(nCycles),
    HHIACT3  = numeric(nCycles),
    HHIECT2  = numeric(nCycles),
    HHIECT3  = numeric(nCycles),
    HHIHPT1  = numeric(nCycles),
    HHIHPT2  = numeric(nCycles),
    HHIHPT3  = numeric(nCycles),
    HHIHPT4  = numeric(nCycles),
    HHISeedlings = numeric(nCycles),
    HHIF1    = numeric(nCycles),
    He_chr6 = numeric(nCycles),
    allelesSI= numeric(nCycles),
    
    # Add dynamic columns for each stage's meanG and varG
    matrix(numeric(nCycles * length(meanG_cols)), ncol = length(meanG_cols), dimnames = list(NULL, meanG_cols)),
    matrix(numeric(nCycles * length(varG_cols)),  ncol = length(varG_cols),  dimnames = list(NULL, varG_cols))
  )
  


  # ---- Create initial parents ----
  source(file = "CreateParents.R")

  # ---- Fill breeding pipeline with unique individuals from initial parents ----
  source(file = "FillPipeline.R")

  # ---- Simulate year effects ----
  P = runif(nCycles)
  

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
    


    for (stage in stages) {
      stage_obj <- get(stage)  # Get the stage object by name
      
      # Calculate mean genetic value
      mean_g <- meanG(stage_obj)
      # Calculate genetic variance
      var_g <- varG(stage_obj)
      
      # Store in output using dynamic column names
      output[year, paste0("meanG_", stage)] <- mean_g
      output[year, paste0("varG_", stage)]  <- var_g
    }
    
    
    
    
    # Inside your loop over years (e.g. for (year in 1:n_years))
    for(stage in stages){
      # Get the population object by name
      pop <- get(stage)
      
      # Get mother and father IDs as characters
      ids <- c(as.character(pop@mother), as.character(pop@father))
      
      # Calculate relative frequencies
      freq <- table(ids) / length(ids)
      
      # Compute Herfindahl–Hirschman Index (HHI)
      hhi_raw <- sum(freq^2)
      
      # Store the HHI value in the output list
      output[[paste0("HHI", stage)]][year] <- hhi_raw
    }
    
    
    pollenHaplo <-  pullMarkerHaplo(Seedlings,
                                    markers = prms$SIPos,
                                    haplo = 2)
    pistilHaplo <-  pullMarkerHaplo(Seedlings,
                                    markers = prms$SIPos,
                                    haplo = 1)
    combinedHaplo <- rbind(pistilHaplo, pollenHaplo)
    # To characters
    haploStrings <- apply(combinedHaplo, 1, paste, collapse = "_")

    # Unique count
    n_unique_haplotypes <- length(unique(haploStrings))
    output$allelesSI[year]=n_unique_haplotypes
    
    
    # Extract genotypes in chr 6
    geno_chr6 <- pullSnpGeno(Seedlings, snpChip = 1, chr = 6, simParam = SP)
    
    # 
    # Allelic frequency per snp
    allele_freqs <- colMeans(geno_chr6) / 2  # Si genotipo codificado 0,1,2 para número de alelos alternativos
    
    # Expected heterocigosity Hs = 2p(1-p)
    Hs <- 2 * allele_freqs * (1 - allele_freqs)
    
    # Mean
    mean_Hs <- mean(Hs, na.rm = TRUE)
    
    #Save
    output$He_chr6[year] = mean_Hs
    

    
  }

  save.image(paste0("../burn_in_folder/Burnin_", REP, ".RData"))
  # Save results from current replicate
  # results = append(results, list(output))
}



# Save results
# saveRDS(results, file = paste0(scenarioName,REP,".rds"))


# ---- Analyze results ----
# source(file = "ANALYZERESULTS.R")


