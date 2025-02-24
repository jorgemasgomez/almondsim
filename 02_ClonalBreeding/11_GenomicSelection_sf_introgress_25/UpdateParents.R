# Update parents

# Replace parents with the lowest ebv, with the new ones cohort


# Update parents

# Replace parents with the lowest ebv, with the new ones cohort


poplist=list(Parents, ECT1)
pop_parents=mergePops(poplist)

#During burning
if (year<=nBurnin){
  Parents = selectInd(pop_parents, nInd = nParents, use = "pheno")
}



#First year after burning 
if (year==nBurnin+1){
  pop_parents = setEBV(pop_parents,gsmodel,value = "bv")
  Parents = selectInd(pop_parents, nInd = nParents, use = "ebv")
  
  pollenHaplo <-  pullMarkerHaplo(Parents,
                                  markers = prms$SIPos,
                                  haplo = "all")
  #We assign randomly to one parent the selfcompatible haplotype
  pollenHaplo[random_selfcompatible_parent:random_selfcompatible_parent,1:length(pollenHaplo[1,])] <- as.matrix(self_compatible_allele)
  
  Parents=setMarkerHaplo(pop = Parents, haplo = pollenHaplo)
  
}

# #Second year and forward after burning

if (year>nBurnin+1){
  pop_parents = setEBV(pop_parents,gsmodel,value = "bv")
  Parents=selectind_selfcompatible(original_pop = pop_parents, haploposition = prms$SIPos,
                                   self_compatible_allele = as.matrix(self_compatible_allele), nindselectedtotal = nParents,
                                   percentage_selfcomp = 0.5, use="ebv" )
  
  
}






