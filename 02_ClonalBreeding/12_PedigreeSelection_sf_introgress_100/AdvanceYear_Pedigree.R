# Advance year

# Advance breeding program by 1 year
# Works backwards through pipeline to avoid copying data



# Stage 12
ECT3 = setPheno(ECT2, varE = VarE, reps = repECT, p = P[year],h2 = c(0.5))

# Stage 11
ECT2 = ECT1

# Stage 10
#TEST_JORGE
ACT3@ebv <- as.matrix(EBV[(length(EBV) - HPT4@nInd - ACT3@nInd + 1):(length(EBV) - HPT4@nInd)])
output$accSel[year] = cor(gv(ACT3), ebv(ACT3))

ECT1= ACT3
ECT1 = selectind_selfcompatible(original_pop = ECT1, haploposition = prms$SIPos,
                                self_compatible_allele = self_compatible_allele,nindselectedtotal = nClonesECT,
                                percentage_selfcomp = 0.5, use="ebv")



# ECT1 = selectInd(ACT3, nInd = nClonesECT, use = "pheno")


# Stage 9
ACT3 = setPheno(ACT2, varE = VarE, reps = repACT, p = P[year], h2 = c(0.4))

# Stage 8
ACT2 = ACT1

# Stage 7
# Use pedigree estimated breeding values to select seedlings for further evaluation
HPT4@ebv= as.matrix(tail(EBV, HPT4@nInd))

output$accSel[year] = cor(gv(HPT4), ebv(HPT4))
ACT1=HPT4
ACT1 = selectind_selfcompatible(original_pop = ACT1, haploposition = prms$SIPos,
                                self_compatible_allele = self_compatible_allele,
                                nindselectedtotal = nClonesACT,
                                percentage_selfcomp = 0.5, use="ebv")


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
F1 = randCrossGamSI_selfcomp_compulsory(prms = prms,pop = Parents, nCrosses = nCrosses, nProgeny = nProgeny,
                                        percentage_selfcompatible = 1)
Parents = setPheno(Parents, varE = VarE, reps = repECT, p = P[year], h2 = c(0.5))

























