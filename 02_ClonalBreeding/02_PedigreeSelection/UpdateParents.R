# Update parents

# Replace parents with the lowest ebv, with the new ones cohort


if (year >= nBurnin+1) {
  poplist=list(Parents, ECT1)
  pop_parents=mergePops(poplist)
  
  pop_parents@ebv <- as.matrix(EBV2$ebv[match(pop_parents@id, EBV2$id)])
  Parents = selectInd(pop_parents, nInd = nParents, use = "ebv")
} else {
  poplist=list(Parents, ECT1)
  pop_parents=mergePops(poplist)
  Parents = selectInd(pop_parents, nInd = nParents, use = "pheno")
} 

