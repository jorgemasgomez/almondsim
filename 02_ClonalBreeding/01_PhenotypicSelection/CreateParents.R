# Create founders

# # Create founder population
# #founderPop = runMacs2(nInd     = nParents,
#                       nChr     = nChr,
#                       segSites = nQtl + nSnp,
#                       genLen   = genLen,
#                       mutRate  = mutRate)
# 
# # Set simulation parameters


# Poner aqui la importacion de los parentales

#La posicion genetica hacerla funcion de la posici?n f?sica (recomb rate del chr  x posici?n fisica)


#Import physical map

genMap_test <- scan("Chr_test_position.txt")
recombination_rate_test<- 1.97e-08
# Convert to genetic map
genMap_test<- genMap_test*recombination_rate_test
genMap=list(genMap_test)

#Import gmatrix for haplotypes
data <- scan("Chr_test_gmatrix.txt")
chr_test <- matrix(data, ncol = 40, byrow = TRUE)
chr_test<-t(chr_test)
haplotypes = list(chr_test)

founderPop = newMapPop(genMap=genMap, haplotypes=haplotypes)

SP = SimParam$new(founderPop)


# Add SNP chip
SP$restrSegSites(nQtl,nSnp)
if (nSnp > 0) {
  SP$addSnpChip(nSnp)
}

# Add traits: trait represents yield
SP$addTraitADG(nQtlPerChr = nQtl,
               mean       = initMeanG,
               var        = initVarG,
               varGxE     = initVarGE)

# Collect pedigree
SP$setTrackPed(TRUE)

# Create founder parents
Parents = newPop(founderPop)

# Set a phenotype to founder parents
Parents = setPheno(Parents, varE = VarE, reps = repECT)
rm(founderPop)
