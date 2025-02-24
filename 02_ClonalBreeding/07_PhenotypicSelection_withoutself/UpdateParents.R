# Update parents

# Replace parents with the lowest ebv, with the new ones cohort


poplist=list(Parents, ECT1)
pop_parents=mergePops(poplist)
Parents = selectInd(pop_parents, nInd = nParents, use = "pheno")
