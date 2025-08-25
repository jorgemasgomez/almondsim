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


# Inicializar listas vacías
genMap <- list()
haplotypes <- list()
# Definir tasa de recombinación como una lista con 8 valores
recombination_rate_test <- c(1.97036e-08, 1.51361e-08, 1.886e-08, 2.22347e-8,2.54212e-8,1.9829e-8,2.61388e-8,2.11591e-8)

# Iterar sobre los cromosomas del 1 al 8
for (i in 1:8) {
  # Leer el archivo de posiciones
  genMap_temp <- scan(paste0("../Haplotypes/Chr_", i, "_position.txt"))
  
  # Convertir a mapa genético utilizando la tasa de recombinación correspondiente
  genMap_temp <- genMap_temp * recombination_rate_test[i]
  
  # Guardar en la lista
  genMap[[i]] <- genMap_temp
  
  # Importar matriz de haplotipos
  data <- scan(paste0("../Haplotypes/Chr_", i, "_gmatrix.txt"))
  chr_temp <- matrix(data, ncol = 60, byrow = TRUE)
  chr_temp <- t(chr_temp)
  
  # Guardar en la lista de haplotipos
  haplotypes[[i]] <- chr_temp
}

# Ahora genMap y haplotypes contienen los datos de los cromosomas 1 a 8


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

#set positions SI locus
#We add +1 for self-compatible allele
nmarkers_sloci <- ceiling(log2(nlocicompatibility+1))


prefix <- sub("_.*", "", locuscompt_position)
base_number <- as.numeric(sub(".*_(\\d+)$", "\\1", locuscompt_position))
secuencia <- base_number + 1:nmarkers_sloci
locuscompt_position_vector <- paste(prefix, secuencia, sep = "_")

## Makes recombination among positions 0
chromosomeSI= prefix
SIPosStart =base_number
SIPosStop=base_number+ nmarkers_sloci

genMap <- SP$genMap
genMap[[chromosomeSI]][SIPosStart:SIPosStop] <-
  genMap[[chromosomeSI]][SIPosStart]
SP$switchGenMap(genMap)

# Create founder parents
Parents = newPop(founderPop, simParam=SP)

#SET SLOCI


slocihaplos<-pullMarkerHaplo(Parents,markers = locuscompt_position_vector, haplo = "all")


# Generate all possible binary combinations
# For each number from 0 to 2^nmarkers_sloci - 1, convert to binary and take only the first nmarkers_sloci bits
combinations <- sapply(0:(2^nmarkers_sloci - 1), function(x) {
  as.numeric(intToBits(x))[1:nmarkers_sloci]
})
# Transpose the combinations and convert them into a data frame
possible_haplotypes <- as.data.frame(t(combinations))
colnames(possible_haplotypes) <- locuscompt_position_vector
# Keep only the first nloci rows

possible_haplotypes <- possible_haplotypes[1:nlocicompatibility, ]

#If you want to introduce self-compatible allele in the population
# possible_haplotypes <- possible_haplotypes[1:nlocicompatibility+1, ]


self_compatible_allele<-as.data.frame(t(combinations))[nlocicompatibility+1,]

#RANDOM ASSIGNMETN

# Extract unique IDs from the row names of slocihaplos
rownames_split <- strsplit(rownames(slocihaplos), "_")  # Split row names by "_"
ids <- sapply(rownames_split, `[`, 1)  # Extract the ID (the part before "_")
unique_ids <- unique(ids)  # Get unique IDs

# Create a copy of slocihaplos to populate with new haplotypes
new_slocihaplos <- slocihaplos

# Loop through each unique ID to assign haplotypes
for (id in unique_ids) {
  # Find rows corresponding to the current ID (e.g., "1_1", "1_2")
  rows_for_id <- which(ids == id)
  
  # Randomly select two different rows from possible_haplotypes
  selected_haplotypes <- possible_haplotypes[sample(1:nrow(possible_haplotypes), 2, replace = FALSE), ]
  
  # Assign the selected haplotypes to the corresponding rows in new_slocihaplos
  new_slocihaplos[rows_for_id[1], ] <- as.numeric(selected_haplotypes[1, ])
  new_slocihaplos[rows_for_id[2], ] <- as.numeric(selected_haplotypes[2, ])
}



Parents=setMarkerHaplo(pop = Parents, haplo = new_slocihaplos)

# Set a phenotype to founder parents
Parents = setPheno(Parents, varE = VarE, reps = repECT)
rm(founderPop)
