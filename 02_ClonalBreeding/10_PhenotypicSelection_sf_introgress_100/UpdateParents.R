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
  Parents = selectInd(pop_parents, nInd = nParents, use = "pheno")
  
  pollenHaplo <-  pullMarkerHaplo(Parents,
                                  markers = prms$SIPos,
                                  haplo = "all")
  #We assign randomly to one parent the selfcompatible haplotype
  pollenHaplo[random_selfcompatible_parent:random_selfcompatible_parent,1:length(pollenHaplo[1,])] <- as.matrix(self_compatible_allele)
  
  Parents=setMarkerHaplo(pop = Parents, haplo = pollenHaplo)
  
}

# #Second year and forward after burning

if (year>nBurnin+1){

  Parents=selectind_selfcompatible(original_pop = pop_parents, haploposition = prms$SIPos,
                           self_compatible_allele = as.matrix(self_compatible_allele), nindselectedtotal = nParents,
                           percentage_selfcomp = 0.5 )
  
  
}






# random_selfcompatible_parent=1
# pollenHaplo <-  pullMarkerHaplo(pop_parents,
#                                 markers = prms$SIPos,
#                                 haplo = "all")
# #We assign randomly to one parent the selfcompatible haplotype
# pollenHaplo[random_selfcompatible_parent:random_selfcompatible_parent,1:length(pollenHaplo[1,])] <- as.matrix(self_compatible_allele)
# 
# pop_parents=setMarkerHaplo(pop = pop_parents, haplo = pollenHaplo)
# 
# 





