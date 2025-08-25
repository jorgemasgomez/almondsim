# Script name: Phenotypic tea clonal breeding program
#
# Authors: Initially developed by Nelson Lubanga, Gregor Gorjanc, Jon Bancic;
# exanded/polished for this publication by Jon Bancic, Philip Greenspoon,
# Chris Gaynor, Gregor Gorjanc
#
# Date Created: 2023-12-06

# ---- Clean environment and load packages ----

# install.packages(pkgs = "AlphaSimR")
library(package = "AlphaSimR")

# pipeline<-FALSE

if (pipeline) {
  scenarioName = "Pheno_withoutself"
  # ---- Create list to store results from reps ----
  results = list()
  
  # ---- Future phase ----
  for(year in (nBurnin+1):(nBurnin+nFuture)) {
    cat("  Working on future year:",year,"\n")
    source(file = "UpdateParents.R")  # Pick parents
    source(file = "AdvanceYear.R")    # Advances yield trials by a year and collects records
    # Report results
    output$meanG[year] = meanG(Seedlings)
    output$varG[year]  = varG(Seedlings)
    # Obtener IDs combinados (madre + padre)
    act1_ids <- c(as.character(ACT1@mother), as.character(ACT1@father))
    ect1_ids <- c(as.character(ECT1@mother), as.character(ECT1@father))
    
    # Obtener IDs combinados (madre + padre)
    act1_ids <- c(as.character(ACT1@mother), as.character(ACT1@father))
    ect1_ids <- c(as.character(ECT1@mother), as.character(ECT1@father))
    # Calcular HHI
    act1_freq <- table(act1_ids) / length(act1_ids)
    act1_hhi_raw <- sum(act1_freq^2)
    act1_hhi_norm <- act1_hhi_raw
    # act1_hhi_norm <- act1_hhi_raw / (1 / (nParents*2))
    
    ect1_freq <- table(ect1_ids) / length(ect1_ids)
    ect1_hhi_raw <- sum(ect1_freq^2)
    ect1_hhi_norm <- ect1_hhi_raw
    # ect1_hhi_norm <- ect1_hhi_raw / (1 / (nParents*2))
    
    
    output$HHIact1[year]=act1_hhi_norm
    output$HHIect1[year]=ect1_hhi_norm
    
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
  
  # Save results from current replicate
  results = append(results, list(output))
  
  
  
  
  
  
} else {
  rm(list = ls())
  source(file = "compatible_crosses.R")
  # ---- Load global parameters ----
  source(file = "GlobalParameters.R")
  scenarioName = "Pheno"
  
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
      source(file = "AdvanceYear.R")    # Advances yield trials by a year and collects records
      # Report results
      output$meanG[year] = meanG(Seedlings)
      output$varG[year]  = varG(Seedlings)
    }
    
    # ---- Future phase ----
    for(year in (nBurnin+1):(nBurnin+nFuture)) {
      cat("  Working on future year:",year,"\n")
      source(file = "UpdateParents.R")  # Pick parents
      source(file = "AdvanceYear.R")    # Advances yield trials by a year and collects records
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








