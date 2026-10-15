# Global Parameters

# ---- Number of simulation replications and breeding cycles ----
nReps = 10 # Independent simulation replicates
nBurnin = 40     # Number of years in burnin phase
nFuture = 40 # Number of years in future phase
startRecords = 35 # Year when training and pedigree record collecting begins
nCycles = nBurnin + nFuture

# ---- Genome simulation ----
nChr    = 8      # Eight chromosomes imported from the external bank
nQtl    = 160    # Per chromosome: 1280 QTL across eight chromosomes
nSnp    = 600    # Per chromosome: 4800 SNP-chip markers across eight chromosomes
genLen  = 1       # Genetic length
PhyLen  = 3e+07   # Physical length
mutRate = 1e-08 # Mutation rate
nlocicompatibility=50
locuscompt_position = "6_5000"
# ---- Initial parents mean and variance ----
initMeanG = 2500    # Phenotypic mean
initVarG  = 150000  # Genetic variance
initVarGE = 150000  # Genotype-by-year interaction variance
VarE      = 2800000 # Single variance

# ---- Breeding program details ----
nParents = 32 # Diploid founders sampled per replicate
nCrosses   = 20  # Number of crosses
nProgeny   = 25   # Number of progenies per cross
nClonesACT = 25  # Number of individuals selected at ACT stage
nClonesECT = 8   # Number of individuals selected at ECT stage
nClonesECTGS = 13
random_selfcompatible_parent=sample(1:nParents*2, 1)
# Effective replication of yield trials
repHPT = 1 # PT aggregate h2 = 0.20
repACT = 1 # AT aggregate h2 = 0.40
repECT = 1 # ET aggregate h2 = 0.50
