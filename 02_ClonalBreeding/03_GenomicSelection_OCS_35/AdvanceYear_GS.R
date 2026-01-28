# Advance year

# Advance breeding program by 1 year
# Works backwards through pipeline to avoid copying data


# Stage 12
ECT3 = setPheno(ECT2, varE = VarE, reps = repECT, p = P[year],h2 = c(0.5))

# Stage 11
ECT2 = ECT1

# Stage 10
ETC1=setEBV(ACT3, gsmodel, value = "bv")
output$accSel[year] = cor(gv(ETC1), ebv(ETC1)) # accuracy based on 2000 inds

### OCS
phen_df <- create_phen_df(ETC1, breed_name = "Almond", sex = NA, herd = NA)
files <- export_snp_by_chromosome(ETC1, out_dir = ".")
sKin <- segIBD(files, map)
files <- unlist(files)
file.remove(files)

cand  <- candes(phen=phen_df, sKin=sKin, cont=cont)
# con <- list(
#   ub.sKin = 1 - (1 - cand$mean$sKin) * (1 - 1/(2*Ne))^(1/L)
# )

con <- list(
  ub.sKin = cand$mean$sKin + cand$mean$sKin*percentage_skin
)



candidates_list <- opticont("max.BV", cand, con, trace=FALSE)


top_list <- candidates_list$parent$Indiv[
  order(candidates_list$parent$oc, decreasing = TRUE)
][1:nClonesECT]

### OCS

ECT1 = selectInd(ETC1, nInd = nClonesECT, candidates = top_list, use = "ebv")

# Stage 9
ACT3 = setPheno(ACT2, varE = VarE, reps = repACT, p = P[year], h2 = c(0.4))

# Stage 8
ACT2 = ACT1

# Stage 7

ACT1 = setEBV(HPT4, gsmodel, value = "bv")
output$accSel[year] = cor(gv(ACT1), ebv(ACT1)) # accuracy based on 2000 inds
ACT1 = selectInd(ACT1, nInd = nClonesACT, use = "ebv")

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

######################################### OCS ##############################

phen_df <- create_phen_df(Parents, breed_name = "Almond", sex = NA, herd = NA)
files <- export_snp_by_chromosome(Parents, out_dir = ".")
sKin <- segIBD(files, map)
files <- unlist(files)
file.remove(files)

cand  <- candes(phen=phen_df, sKin=sKin, cont=cont)
# con <- list(
#   ub.sKin = 1 - (1 - cand$mean$sKin) * (1 - 1/(2*Ne))^(1/L)
# )

con <- list(
  ub.sKin = cand$mean$sKin + cand$mean$sKin*percentage_skin
)
candidates_list <- opticont("max.BV", cand, con, trace=FALSE)

Candidate <- candidates_list$parent[,  c("Indiv", "Sex", "oc", "herd")]
Candidate[Candidate$oc>0.001,]

Candidate$Sex <- "male"
Candidate$n <- noffspring(Candidate, N=nProgeny*nCrosses/2)$nOff
Candidate$Sex <- NA
Mating <- matings_hermaf(Candidate, Kin=sKin,pop = Parents, SIPos = prms$SIPos)

Parent1 <- pmin(Mating$Sire, Mating$Dam)
Parent2 <- pmax(Mating$Sire, Mating$Dam)

# Create df with pairs
Mating_ord <- data.frame(Parent1, Parent2, n = Mating$n)

# accumulate pairs
Mating_acum <- aggregate(n ~ Parent1 + Parent2, data = Mating_ord, sum)

impossible_crosses <- getImpossibleCrosses(Parents, SIPos= prms$SIPos)

impossible_crosses$Cross <- "Impossible"


# Merge
Mating_acum <- merge(
  Mating_acum,
  impossible_crosses[, c("Parent1","Parent2","Cross")],
  by = c("Parent1","Parent2"),
  all.x = TRUE
)

######################################### OCS ##############################

prms<- list(SIPos=locuscompt_position_vector, random_failure=FALSE)

F1 <- NULL

for (i in seq_len(nrow(Mating_acum))) {
  
  if (is.na(Mating_acum$n[i]) || Mating_acum$n[i] <= 0) next
  
  cross_plan <- matrix(
    c(Mating_acum$Parent1[i],
      Mating_acum$Parent2[i]),
    ncol = 2,
    byrow = TRUE
  )
  
  F1_i <- randCrossGamSI(
    prms       = prms,
    pop        = Parents,
    cross_plan = cross_plan,
    nProgeny   = Mating_acum$n[i],
    nCrosses = 1
  )
  
  F1 <- if (is.null(F1)) F1_i else c(F1, F1_i)
}


Parents = setPheno(Parents, varE = VarE, reps = repECT, p = P[year], h2 = c(0.5))




