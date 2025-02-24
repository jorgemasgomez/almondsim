# Advance year

# Advance breeding program by 1 year
# Works backwards through pipeline to avoid copying data


# Stage 12
ECT3 = setPheno(ECT2, varE = VarE, reps = repECT, p = P[year],h2 = c(0.5))

# Stage 11
ECT2 = ECT1

# Stage 10
output$accSel[year] = cor(gv(ACT3), pheno(ACT3)) # accuracy based on 2000 inds
ECT1 = selectInd(ACT3, nInd = nClonesECT, use = "pheno")

# Stage 9
ACT3 = setPheno(ACT2, varE = VarE, reps = repACT, p = P[year], h2 = c(0.4))

# Stage 8
ACT2 = ACT1

# Stage 7
output$accSel[year] = cor(gv(HPT4), pheno(HPT4)) # accuracy based on 2000 inds
ACT1 = selectInd(HPT4, nInd = nClonesACT, use = "pheno")

# Stage 6
HPT4 = setPheno(HPT3, varE = VarE, reps = repHPT, p = P[year], h2 = c(0.2))

# Stage 5
HPT3 = HPT2

# Stage 4
HPT2 = HPT1

# Stage 3
HPT1 = Seedlings

# Stage 2
Seedlings = F1

# Stage 1
# Crossing block
# F1 = randCross(Parents, nCrosses = nCrosses, nProgeny = nProgeny)
prms<- list(SIPos=locuscompt_position_vector, random_failure=FALSE)
F1 = randCrossGamSI(prms = prms,pop = Parents, nCrosses = nCrosses, nProgeny = nProgeny)
