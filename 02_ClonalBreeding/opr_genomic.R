# Sequential OPR: test candidates in descending GEBV order.
# Reference diversity is fixed at the start of EACH year, not each swap.
# All SNPs, including monomorphic markers, are kept in the metric denominator.
opr_expected_heterozygosity <- function(dosage_sum,n) {
  frequency <- dosage_sum/(2*n)
  mean(2*frequency*(1-frequency))
}

opr_sequential <- function(geno,bv,old,max_new,max_loss,max_candidates=50L) {
  stopifnot(nrow(geno)==length(bv),ncol(geno)>0L,all(is.finite(bv)),
            !anyNA(geno),all(geno>=0 & geno<=2),!anyDuplicated(old),
            max_new>=0L,max_new<=length(old),max_loss>=0,max_loss<=1,max_candidates>=0L)
  n <- length(old); chosen <- old
  sums <- colSums(geno[chosen,,drop=FALSE])
  before <- opr_expected_heterozygosity(sums,n)
  minimum <- (1-max_loss)*before
  candidate_ids <- setdiff(seq_along(bv),old)
  candidate_ids <- candidate_ids[order(-bv[candidate_ids],candidate_ids)]
  candidate_ids <- head(candidate_ids,max_candidates)
  trace <- list();new <- 0L;rejected <- 0L
  for(candidate in candidate_ids) {
    if(new>=max_new) break
    # Exact equivalent of merge parents + candidate, retain best n by GEBV.
    # A tie favors the current parent, so a no-op never counts as replacement.
    worst <- which.min(bv[chosen]);removed <- chosen[worst]
    if(bv[candidate]<=bv[removed]) {
      trace[[length(trace)+1L]] <- data.frame(candidate=candidate,removed=NA_integer_,
                                            accepted=FALSE,reason='not_better',He=NA_real_)
      next
    }
    proposed_sums <- sums-geno[removed,]+geno[candidate,]
    proposedHe <- opr_expected_heterozygosity(proposed_sums,n)
    accept <- proposedHe>=minimum-1e-12
    trace[[length(trace)+1L]] <- data.frame(candidate=candidate,removed=removed,
                                          accepted=accept,reason=if(accept) 'accepted' else 'diversity',He=proposedHe)
    if(accept) {
      chosen[worst] <- candidate;sums <- proposed_sums
      new <- sum(!chosen %in% old)
    } else rejected <- rejected+1L
  }
  after <- opr_expected_heterozygosity(sums,n)
  stopifnot(new<=max_new,after>=minimum-1e-12,length(unique(chosen))==n)
  list(ids=chosen,new=new,beforeHe=before,afterHe=after,minHe=minimum,
       lossPercent=if(before>0) 100*(before-after)/before else 0,
       expectedBV=mean(bv[chosen]),tested=length(trace),rejectedDiversity=rejected,
       trace=if(length(trace)) do.call(rbind,trace) else data.frame())
}

redesign_opr <- function(state,candidates,model) {
  pool <- redesign_predict(redesign_merge_latest(list(state$parents,candidates)),'GS',model)
  old <- match(as.character(state$parents@id),as.character(pool@id))
  result <- opr_sequential(pullSnpGeno(pool,simParam=SP),as.numeric(ebv(pool)),old,
                           length(old),state$config$opr_max_diversity_loss,
                           state$config$opr_candidate_limit)
  result$parents <- pool[result$ids]
  stopifnot(nInd(result$parents)==state$config$parents,
            result$new==sum(!result$parents@id %in% state$parents@id))
  result
}
