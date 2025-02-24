# Global Parameters

# ---- Number of simulation replications and breeding cycles ----
nReps   = 1    # Number of simulation replicates
nBurnin = 40      # Number of years in burnin phase
nFuture = 40      # Number of years in future phase
startRecords = 35 # Year when training and pedigree record collecting begins
nCycles = nBurnin + nFuture

# ---- Genome simulation ----

nQtl    = 160     # number of qtl per chromosome: 15 chr x 160 qtl = 2400 qtls
nSnp    = 600     # simulate snp chip with 9000 markers
genlen  = 1       # genetic length
phylen  = 1e+08   # physical length
mutrate = 2.5e-08 # mutation rate
nlocicompatibility=50
locuscompt_position = "1_5000"

# ---- Initial parents mean and variance ----
initMeanG = 2500    # Phenotypic mean
initVarG  = 150000  # Genetic variance
initVarGE = 150000  # Genotype-by-year interaction variance
VarE      = 2800000 # Single variance

# ---- Breeding program details ----
nParents   = 30   # Number of parents (and founders)
nCrosses   = 20  # Number of crosses
nProgeny   = 50   # Number of progenies per cross
nClonesACT = 50  # Number of individuals selected at ACT stage
nClonesECT = 3   # Number of individuals selected at ECT stage

# Effective replication of yield trials
repHPT     = 1   # h2 = 0.2
repACT     = 1  # h2 = 0.4
repECT     = 1  # h2 = 0.5