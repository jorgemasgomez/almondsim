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
  scenarioName = "GS_OCS_25"
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
    # Convertir a caracteres para comparar fácilmente
    haploStrings <- apply(combinedHaplo, 1, paste, collapse = "_")
    
    # Contar únicos
    n_unique_haplotypes <- length(unique(haploStrings))
    output$allelesSI[year]=n_unique_haplotypes
    
    # Extraer genotipos SNP del cromosoma 6
    geno_chr6 <- pullSnpGeno(Seedlings, snpChip = 1, chr = 6, simParam = SP)
    
    # La matriz geno_chr6 tiene filas = individuos, columnas = SNPs
    # Calcular frecuencia alélica para cada SNP
    allele_freqs <- colMeans(geno_chr6) / 2  # Si genotipo codificado 0,1,2 para número de alelos alternativos
    
    # Heterocigosidad esperada Hs = 2p(1-p)
    Hs <- 2 * allele_freqs * (1 - allele_freqs)
    
    # Resumen de la diversidad en cromosoma 6
    mean_Hs <- mean(Hs, na.rm = TRUE)
    
    # Mostrar resultado
    output$He_chr6[year] = mean_Hs
    
    }
  
  
    # Guarda resultados del año actual
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





