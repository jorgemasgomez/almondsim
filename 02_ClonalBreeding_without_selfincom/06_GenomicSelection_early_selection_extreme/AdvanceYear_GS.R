# Advance year

# Advance breeding program by 1 year
# Works backwards through pipeline to avoid copying data



# Stage 9
ACT3 = setPheno(ACT2, varE = VarE, reps = repACT, p = P[year], h2 = c(0.4))

# Stage 8
ACT2 = ACT1

# Stage 7

ACT1 = setEBV(HPT4, gsmodel, value = "bv")
# output$accSel[year] = cor(gv(ACT1), ebv(ACT1)) # accuracy based on 2000 inds
ACT1 = selectInd(ACT1, nInd = nClonesACT, use = "ebv")

# Stage 6
HPT4 = setPheno(HPT3, varE = VarE, reps = repHPT, p = P[year], h2 = c(0.2))

# Stage 5
HPT3 = HPT2

# Stage 4
HPT2 = HPT1

# Stage 3
HPT1 = setEBV(Seedlings, gsmodel, value = "bv")
output$accSel[year] = cor(gv(HPT1), ebv(HPT1)) # accuracy based on 2000 inds
HPT1 = selectInd(HPT1, nInd = nClonesHPT, use = "ebv")

# Stage 2
Seedlings = F1

# Stage 1
# Crossing block
F1 = randCross(Parents, nCrosses = nCrosses, nProgeny = nProgeny)

# prms<- list(SIPos=locuscompt_position_vector, random_failure=FALSE)
# F1 = randCrossGamSI(prms = prms,pop = Parents, nCrosses = nCrosses, nProgeny = nProgeny)
Parents = setPheno(Parents, varE = VarE, reps = repECT, p = P[year], h2 = c(0.5))




