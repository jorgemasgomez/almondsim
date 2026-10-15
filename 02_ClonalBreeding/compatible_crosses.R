# Mother/father pairs remain fixed while incompatible pollen is rejected.
# Each crossing event produces its full, equal-sized family or fails explicitly.
randCrossGamSI <- function(prms,pop,nCrosses=NA,nProgeny,cross_plan=NA) {
 state <- getOption('almond.cross_state')
 allowed <- ocs_compatibility(pop,prms$SIPos,self_compatible_allele)
 if(!is.null(state)) {
  age <- state$year-state$birth[as.character(pop@id)]
  if(anyNA(age)) stop('Unknown parental birth year in the crossing plan.')
  allowed <- allowed & outer(age>=state$config$minimum_mother_age,
                             age>=state$config$minimum_pollen_age,'&')
 }
 diag(allowed) <- FALSE
 if(all(is.na(cross_plan))) {
  if(length(nCrosses)!=1L||is.na(nCrosses)||nCrosses<1L) stop('A positive crossing budget is required.')
  pairs <- which(allowed,arr.ind=TRUE)
  if(!nrow(pairs)) stop('No S-compatible, non-self pair satisfies the parental age limits.')
  plan <- pairs[sample.int(nrow(pairs),nCrosses,replace=TRUE),,drop=FALSE]
 } else {
  plan <- as.matrix(cross_plan)
  if(is.character(plan)) plan <- matrix(match(plan,as.character(pop@id)),ncol=2)
  if(ncol(plan)!=2L||anyNA(plan)||any(plan<1L)||any(plan>nInd(pop))||!all(allowed[plan])) stop('Invalid S-locus or age-constrained crossing plan.')
  if(!is.na(nCrosses)&&nrow(plan)!=nCrosses) stop('Crossing plan does not match the budget.')
 }
 maternal1 <- pullMarkerHaplo(pop,markers=prms$SIPos,haplo=1,simParam=SP)
 maternal2 <- pullMarkerHaplo(pop,markers=prms$SIPos,haplo=2,simParam=SP)
 self <- as.numeric(self_compatible_allele)
 families <- lapply(seq_len(nrow(plan)),function(i) {
  kept <- list();remaining <- as.integer(nProgeny)
  for(attempt in seq_len(100L)) {
   if(!remaining) break
   # Never resample parents during rejection: realized contributions stay fixed.
   batch <- makeCross(pop,crossPlan=plan[i,,drop=FALSE],nProgeny=max(4L,2L*remaining),simParam=SP)
   pollen <- pullMarkerHaplo(batch,markers=prms$SIPos,haplo=2,simParam=SP)
   m <- plan[i,1L]
   accepted <- apply(pollen,1,function(h) all(h==self)||
                    (!all(h==maternal1[m,])&&!all(h==maternal2[m,])))
   eligible <- which(accepted)
   if(!length(eligible)) next
   take <- min(remaining,length(eligible))
   kept[[length(kept)+1L]] <- batch[eligible[seq_len(take)]]
   remaining <- remaining-take
  }
  if(remaining) stop('Compatible family incomplete after 100 attempts: ',
    as.character(pop@id[plan[i,1]]),' x ',as.character(pop@id[plan[i,2]]))
  family <- mergePops(kept)
  stopifnot(nInd(family)==nProgeny,all(family@mother==pop@id[plan[i,1]]),all(family@father==pop@id[plan[i,2]]))
  family
 })
 result <- mergePops(families)
 if(!is.null(state)) state$speed_cross_trace <- list(
  mother_ages=state$year-state$birth[as.character(result@mother)],
  father_ages=state$year-state$birth[as.character(result@father)],
  cross_ids=matrix(as.character(pop@id[plan]),ncol=2),cross_count=nrow(plan))
 result
}
