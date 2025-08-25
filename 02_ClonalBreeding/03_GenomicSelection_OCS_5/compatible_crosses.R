##**************************************************************
##******** Function for implementing Self-Incompatibility ******
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
    
    diff1 <- pistilHaplo1 - pollenHaplo
    diff1Mismatch <- apply(X = diff1,
                           FUN = function(x)
                             any(x != 0),MARGIN = 1)
    
    diff2 <- pistilHaplo2 - pollenHaplo
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