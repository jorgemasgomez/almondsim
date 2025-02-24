# ##**************************************************************
# ##******** Function for implementing Self-Incompatibility ******
# ##**************************************************************
# 
# ## Arguments:
# 
# ## prms - Parameter values
# ## pop - the population being crossed
# ## nCrosses - the number of crosses performed
# ## nProgeny - the number of progeny per cross
# ## cross_plan - if applicable, a supplied cross plan for initial crossing
# 
# ## Output:
# 
# ## sward_sel - a list containing the selected swards, updated costs
# ## (if applicable), accuracy measures
# 
# ## Randomly cross within a population with self-incompality alleles
# ## This function assumes that all individuals produce both female and male
# ## gametes (i.e., are dioecious), which is true for red clover.
# 
# randCrossGamSI <- function(prms, pop, nCrosses = NA,
#                            nProgeny,cross_plan = NA) {
#   whileloopcount <- 0
#   while(whileloopcount == 0 || (exists("newPop") && newPop@nInd < nProgeny * nCrosses)){
#     ## These are the crosses that will be attempted, chosen at random
#     
#     if(is.na(cross_plan[1]))
#       candidateOffspring <-
#         randCross(pop,nCrosses = nCrosses,
#                   nProgeny = nProgeny)
#     
#     if(!is.na(cross_plan[1]))
#       candidateOffspring <-
#         makeCross(pop = pop,
#                   crossPlan = cross_plan,
#                   nProgeny = nProgeny)
#     
#     ## The mothers of the candidate offspring
#     mothers <- candidateOffspring@mother
#     
#     ## Get the relevant genotypic information at the SI loci
#     pollenHaplo <-  pullMarkerHaplo(candidateOffspring,
#                                     markers = prms$SIPos,
#                                     haplo = 2)
#     pistilHaplo1 <- pullMarkerHaplo(pop,
#                                     markers = prms$SIPos,
#                                     haplo = 1)
#     pistilHaplo2 <- pullMarkerHaplo(pop,
#                                     markers = prms$SIPos,
#                                     haplo = 2)
#     pistilHaplo1 <- pistilHaplo1[paste(mothers,"1", sep = "_"),]
#     pistilHaplo2 <- pistilHaplo2[paste(mothers,"2", sep = "_"),]
#     
#     ## Use the SI genotypic loci to determine which crosses are
#     ## compatible. The compatible crosses are those for which there is mis-match
#     ## between the pollen haplotype and both pistil haplotypes.
#     
#     diff1 <- pistilHaplo1 - pollenHaplo
#     diff1Mismatch <- apply(X = diff1,
#                            FUN = function(x)
#                              any(x != 0),MARGIN = 1)
#     
#     diff2 <- pistilHaplo2 - pollenHaplo
#     diff2Mismatch <- apply(X = diff2,
#                            FUN = function(x)
#                              any(x != 0),MARGIN =1)
#     compatibleCrosses <- which(diff1Mismatch & diff2Mismatch)
#     
#     
#     
#     ## A test with random failure of crosses
#     if(prms$random_failure)
#       newPop <-
#       selectInd(candidateOffspring,
#                 use = "rand",
#                 nInd = length(compatibleCrosses))
#     
#     ## Select all the offspring resulting from compatible
#     ## crosses
#     if(!prms$random_failure){
#       if(packageVersion("AlphaSimR")<"1.5.3")
#         ## After updating to AlphaSimR version 1.5.3 this no longer is working...
#         newPop_prev <- selectInd(candidateOffspring,
#                             use = "rand",
#                             candidates =
#                               as.character(compatibleCrosses),
#                             nInd = length(compatibleCrosses))
#       
#       if(packageVersion("AlphaSimR")>="1.5.3")
#         ## After updating to AlphaSimR version 1.5.3 this works
#         newPop_prev <- selectInd(candidateOffspring,
#                             use = "rand",
#                             candidates = candidateOffspring@id[compatibleCrosses],
#                             nInd = length(compatibleCrosses))
#         
#     if(whileloopcount==0){
#       newPop = newPop_prev
#       
#     }
#       else{
#         popList = list(newPop, newPop_prev)
#         newPop = mergePops(popList)
#       }
# 
#     
#     
#     whileloopcount=whileloopcount+1
#     
#     if (whileloopcount>10){
#       print("More than 10 iterations to complete compatible progeny")
#       break
#     }
#     }
#     
#   }
#   if(newPop@nInd == 0){
#     return(NA)
#   }else{
#     return(newPop)
#   }
# }

##**************************************************************
##*
##**************************************************************
##******** Function for implementing Self-Incompatibility with one self-compatible allele******
##**************************************************************

## Arguments:

## prms - Parameter values
## pop - the population being crossed
## nCrosses - the number of crosses performed
## nProgeny - the number of progeny per cross
## cross_plan - if applicable, a supplied cross plan for initial crossing

## Output:

## sward_sel - a list containing the selected swards, updated costs
## (if applicable), accuracy measures

## Randomly cross within a population with self-incompality alleles
## This function assumes that all individuals produce both female and male
## gametes (i.e., are dioecious), which is true for red clover.

randCrossGamSI <- function(prms, pop, nCrosses = NA,
                           nProgeny,cross_plan = NA) {
  whileloopcount <- 0
  while(whileloopcount == 0 || (exists("newPop") && newPop@nInd < nProgeny * nCrosses)){
    ## These are the crosses that will be attempted, chosen at random
    
    if(is.na(cross_plan[1]))
      candidateOffspring <-
        randCross(pop,nCrosses = nCrosses,
                  nProgeny = nProgeny)
    
    if(!is.na(cross_plan[1]))
      candidateOffspring <-
        makeCross(pop = pop,
                  crossPlan = cross_plan,
                  nProgeny = nProgeny)
    
    ## The mothers of the candidate offspring
    mothers <- candidateOffspring@mother
    
    ## Get the relevant genotypic information at the SI loci
    pollenHaplo <-  pullMarkerHaplo(candidateOffspring,
                                    markers = prms$SIPos,
                                    haplo = 2)
    pistilHaplo1 <- pullMarkerHaplo(pop,
                                    markers = prms$SIPos,
                                    haplo = 1)
    pistilHaplo2 <- pullMarkerHaplo(pop,
                                    markers = prms$SIPos,
                                    haplo = 2)
    pistilHaplo1 <- pistilHaplo1[paste(mothers,"1", sep = "_"),]
    pistilHaplo2 <- pistilHaplo2[paste(mothers,"2", sep = "_"),]
    
    ## Use the SI genotypic loci to determine which crosses are
    ## compatible. The compatible crosses are those for which there is mis-match
    ## between the pollen haplotype and both pistil haplotypes.
    
    
    self_compatible_allele_matrix <- as.matrix(self_compatible_allele)
    is_ref_match <- apply(pollenHaplo, 1, function(row) all(row == self_compatible_allele))
    
    
    diff1 <- pistilHaplo1 - pollenHaplo
    diff2 <- pistilHaplo2 - pollenHaplo

    #If selfcompatible reference allele in the pollen haplotype we place -2 to avoid the incompatibility
    
    diff1[is_ref_match, ] <- -2
    diff2[is_ref_match, ] <- -2
    
  
    diff1Mismatch <- apply(X = diff1,
                           FUN = function(x)
                             any(x != 0),MARGIN = 1)
    diff2Mismatch <- apply(X = diff2,
                           FUN = function(x)
                             any(x != 0),MARGIN =1)
    
    
    compatibleCrosses <- which(diff1Mismatch & diff2Mismatch)
    
    
    
    ## A test with random failure of crosses
    if(prms$random_failure)
      newPop <-
      selectInd(candidateOffspring,
                use = "rand",
                nInd = length(compatibleCrosses))
    
    ## Select all the offspring resulting from compatible
    ## crosses
    if(!prms$random_failure){
      if(packageVersion("AlphaSimR")<"1.5.3")
        ## After updating to AlphaSimR version 1.5.3 this no longer is working...
        newPop_prev <- selectInd(candidateOffspring,
                                 use = "rand",
                                 candidates =
                                   as.character(compatibleCrosses),
                                 nInd = length(compatibleCrosses))
      
      if(packageVersion("AlphaSimR")>="1.5.3")
        ## After updating to AlphaSimR version 1.5.3 this works
        newPop_prev <- selectInd(candidateOffspring,
                                 use = "rand",
                                 candidates = candidateOffspring@id[compatibleCrosses],
                                 nInd = length(compatibleCrosses))
      
      
      if(whileloopcount==0){
        newPop = newPop_prev
        
      }
      else{
        popList = list(newPop, newPop_prev)
        newPop = mergePops(popList)
      }
      
      
      
      whileloopcount=whileloopcount+1
      
      if (whileloopcount>10){
        print("More than 10 iterations to complete compatible progeny")
        break
      }
    }
    
  }
  if(newPop@nInd == 0){
    return(NA)
  }else{
    return(newPop)
  }
}

##**************************************************************
##*
##**************************************************************
##******** Function for implementing Self-Incompatibility with one self-compatible allele,
##* indicating a percentage compulsory of selfcompatibles seedlings **********************
##**************************************************************

## Arguments:

## prms - Parameter values
## pop - the population being crossed
## nCrosses - the number of crosses performed
## nProgeny - the number of progeny per cross
## cross_plan - if applicable, a supplied cross plan for initial crossing

## Output:

## sward_sel - a list containing the selected swards, updated costs
## (if applicable), accuracy measures

## Randomly cross within a population with self-incompality alleles
## This function assumes that all individuals produce both female and male
## gametes (i.e., are dioecious), which is true for red clover.

randCrossGamSI_selfcomp_compulsory <- function(prms, pop, nCrosses = NA,
                           nProgeny,cross_plan = NA, percentage_selfcompatible=0.25) {
  whileloopcount <- 0
  self_compatible_seedling_number <- 0
  
  while(whileloopcount == 0 || (exists("newPop") && newPop@nInd < nProgeny * nCrosses)){
    ## These are the crosses that will be attempted, chosen at random
    
    if(is.na(cross_plan[1]))
      candidateOffspring <-
        randCross(pop,nCrosses = nCrosses,
                  nProgeny = nProgeny)
    
    if(!is.na(cross_plan[1]))
      candidateOffspring <-
        makeCross(pop = pop,
                  crossPlan = cross_plan,
                  nProgeny = nProgeny)
    
    ## The mothers of the candidate offspring
    mothers <- candidateOffspring@mother
    
    ## Get the relevant genotypic information at the SI loci
    pollenHaplo <-  pullMarkerHaplo(candidateOffspring,
                                    markers = prms$SIPos,
                                    haplo = 2)
    pistilHaplo1 <- pullMarkerHaplo(pop,
                                    markers = prms$SIPos,
                                    haplo = 1)
    pistilHaplo2 <- pullMarkerHaplo(pop,
                                    markers = prms$SIPos,
                                    haplo = 2)
    pistilHaplo1 <- pistilHaplo1[paste(mothers,"1", sep = "_"),]
    pistilHaplo2 <- pistilHaplo2[paste(mothers,"2", sep = "_"),]
    
    ## Use the SI genotypic loci to determine which crosses are
    ## compatible. The compatible crosses are those for which there is mis-match
    ## between the pollen haplotype and both pistil haplotypes.
    
    
    self_compatible_allele_matrix <- as.matrix(self_compatible_allele)
    is_ref_match <- apply(pollenHaplo, 1, function(row) all(row == self_compatible_allele))
    
    
    diff1 <- pistilHaplo1 - pollenHaplo
    diff2 <- pistilHaplo2 - pollenHaplo
    
    #If selfcompatible reference allele in the pollen haplotype we place -2 to avoid the incompatibility
    
    diff1[is_ref_match, ] <- -2
    diff2[is_ref_match, ] <- -2
    
    
    diff1Mismatch <- apply(X = diff1,
                           FUN = function(x)
                             any(x != 0),MARGIN = 1)
    diff2Mismatch <- apply(X = diff2,
                           FUN = function(x)
                             any(x != 0),MARGIN =1)
    
    
    compatibleCrosses <- which(diff1Mismatch & diff2Mismatch)
    
    
    
    ## A test with random failure of crosses
    if(prms$random_failure)
      newPop <-
      selectInd(candidateOffspring,
                use = "rand",
                nInd = length(compatibleCrosses))
    
    ## Select all the offspring resulting from compatible
    ## crosses
    if(!prms$random_failure){
      if(packageVersion("AlphaSimR")<"1.5.3")
        ## After updating to AlphaSimR version 1.5.3 this no longer is working...
        newPop_prev <- selectInd(candidateOffspring,
                                 use = "rand",
                                 candidates =
                                   as.character(compatibleCrosses),
                                 nInd = length(compatibleCrosses))
      
      if(packageVersion("AlphaSimR")>="1.5.3")
        ## After updating to AlphaSimR version 1.5.3 this works
        newPop_prev <- selectInd(candidateOffspring,
                                 use = "rand",
                                 candidates = candidateOffspring@id[compatibleCrosses],
                                 nInd = length(compatibleCrosses))
      
      
      #We prioritize selfcompatibles up to a percentage
      
      if(self_compatible_seedling_number < nCrosses*nProgeny*percentage_selfcompatible){
        
        slocusHaplo <-  pullMarkerHaplo(newPop_prev,
                                        markers = prms$SIPos,
                                        haplo = "all")
        is_ref_match <- apply(slocusHaplo, 1, function(row) all(row == self_compatible_allele))
        
        self_compatible_ind <- which(is_ref_match)
        ids_selfcompatibles_ind<- gsub("(_1|_2)$", "", names(self_compatible_ind))
        ids_selfcompatibles_ind<- unique(ids_selfcompatibles_ind)
        
        newPop_prev=selectInd(newPop_prev, nInd = length(ids_selfcompatibles_ind),
                               candidates = ids_selfcompatibles_ind, use = "rand")
        
        self_compatible_seedling_number=self_compatible_seedling_number + length(newPop_prev@id)
        
      }
      
      
      
      
      
      
      if(whileloopcount==0){
        newPop = newPop_prev
        
      }
      else{
        popList = list(newPop, newPop_prev)
        newPop = mergePops(popList)
      }
      
      
      
      whileloopcount=whileloopcount+1
      
      if (whileloopcount>100){
        print("More than 10 iterations to complete compatible progeny")
        break
      }
    }
    
  }
  if(newPop@nInd == 0){
    return(NA)
  }else{
    return(newPop)
  }
}

##**************************************************************
##*
##*
##*
##*
##*
##*
##*
##*
##*
##*
##*
##*
##*

selectind_selfcompatible<- function(original_pop,haploposition, self_compatible_allele, nindselectedtotal,
                                    percentage_selfcomp=0.5, use="pheno"){
  
  #Check the alleles, introduce those self-compatible up to 50% or if there are less 
  slocusHaplo <-  pullMarkerHaplo(original_pop,
                                  markers = haploposition,
                                  haplo = "all")
  is_ref_match <- apply(slocusHaplo, 1, function(row) all(row == self_compatible_allele))
  
  self_compatible_parents <- which(is_ref_match)
  ids_selfcompatibles_parents <- gsub("(_1|_2)$", "", names(self_compatible_parents))
  ids_selfcompatibles_parents<- unique(ids_selfcompatibles_parents)
  
  
  if (length(ids_selfcompatibles_parents) <= as.integer(nindselectedtotal * percentage_selfcomp)) {
    
    pop1_parents=selectInd(original_pop, nInd = length(ids_selfcompatibles_parents),
                           candidates = ids_selfcompatibles_parents)
    
    
    rest <- setdiff(original_pop@id, pop1_parents@id)
    
    
    pop2_parents=selectInd(original_pop, nInd = nindselectedtotal-length(pop1_parents@id),
                           candidates = rest, use = use)
    
    
    poplist=list(pop1_parents, pop2_parents)
    final_pop=mergePops(poplist)
    
    return(final_pop)
  } else {
    
    pop1_parents=selectInd(original_pop, nInd = as.integer(nindselectedtotal * percentage_selfcomp),
                           candidates = ids_selfcompatibles_parents, use = use )
    
    rest <- setdiff(original_pop@id, pop1_parents@id)
    
    pop2_parents=selectInd(original_pop, nInd = nindselectedtotal-length(pop1_parents@id),
                          candidates = rest, use = use)
    
    poplist=list(pop1_parents, pop2_parents)
    final_pop=mergePops(poplist)
    
    return(final_pop)
  }
  
  
  
}       
