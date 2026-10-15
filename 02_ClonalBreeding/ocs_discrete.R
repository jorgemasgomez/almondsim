# Genomic coancestry, compatibility, metrics and capped local refinement;
# no claim of a globally optimal cardinality-constrained solution.
# Equal family sizes give each parental appearance contribution 1/(2*nCrosses).

ocs_fixed_frequency <- function(pop) {
  geno <- pullSnpGeno(pop,simParam=SP)
  frequency <- colMeans(geno)/2
  if(sum(2*frequency*(1-frequency))<=0) stop('OCS needs polymorphic SNPs.')
  frequency
}

ocs_kinship <- function(pop,frequency) {
  geno <- pullSnpGeno(pop,simParam=SP)
  stopifnot(ncol(geno)==length(frequency),!anyNA(geno))
  centered <- sweep(geno,2,2*frequency,'-')
  K <- tcrossprod(centered)/(2*sum(2*frequency*(1-frequency)))
  dimnames(K) <- list(as.character(pop@id),as.character(pop@id))
  K
}

ocs_compatibility <- function(pop,markers,self_compatible=NULL) {
  h1 <- pullMarkerHaplo(pop,markers=markers,haplo=1,simParam=SP)
  h2 <- pullMarkerHaplo(pop,markers=markers,haplo=2,simParam=SP)
  a <- apply(h1,1,paste,collapse='_');b <- apply(h2,1,paste,collapse='_')
  self <- if(is.null(self_compatible)) NA_character_ else paste(as.numeric(self_compatible),collapse='_')
  # Rows = mother, columns = father. At least one pollen allele must be accepted.
  accepts_a <- outer(a,a,'!=') & outer(b,a,'!=')
  accepts_b <- outer(a,b,'!=') & outer(b,b,'!=')
  if(!is.na(self)) {
    accepts_a[,a==self] <- TRUE
    accepts_b[,b==self] <- TRUE
  }
  allowed <- accepts_a | accepts_b
  diag(allowed) <- FALSE
  dimnames(allowed) <- list(as.character(pop@id),as.character(pop@id))
  allowed
}

ocs_metrics <- function(plan,K,bv) {
  counts <- tabulate(as.integer(plan),nbins=nrow(K))
  contribution <- counts/length(plan)
  active <- which(counts>0)
  list(kinship=as.numeric(crossprod(contribution[active],K[active,active,drop=FALSE]%*%contribution[active])),
       gain=sum(contribution*bv),parents=length(active),counts=counts,
       contribution=contribution,mate_kinship=mean(K[plan]))
}

ocs_improve <- function(plan,K,bv,allowed,max_parents,mode,bound=Inf,sweeps=12L) {
  n <- nrow(K);slots <- length(plan);crosses <- nrow(plan)
  counts <- tabulate(plan,nbins=n)
  for(sweep in seq_len(sweeps)) {
    changes <- 0L
    for(slot in sample.int(slots)) {
      row <- (slot-1L)%%crosses+1L;column <- (slot-1L)%/%crosses+1L
      old <- plan[row,column];mate <- plan[row,3L-column]
      counts[old] <- counts[old]-1L
      active <- which(counts>0)
      options <- which(if(column==1L) allowed[,mate] else allowed[mate,])
      options <- options[counts[options]<ocs_cap_slots(slots/2)]
      if(length(active)>=max_parents) options <- options[counts[options]>0L]
      if(!length(options)) stop('OCS lost a feasible parental assignment.')
      weight <- counts[active]/slots
      Kc <- as.numeric(K[,active,drop=FALSE]%*%weight)
      base <- sum(weight*Kc[active])
      kinship <- base+2*Kc[options]/slots+diag(K)[options]/slots^2
      if(mode=='gain') {
        feasible <- kinship<=bound+1e-10
        options <- options[feasible];kinship <- kinship[feasible]
        if(!length(options)) stop('OCS repair failed to preserve the realized bound.')
        order <- order(-bv[options],kinship,K[options,mate],options)
      } else {
        order <- order(kinship,-bv[options],K[options,mate],options)
      }
      chosen <- options[order[1L]]
      plan[row,column] <- chosen
      counts[chosen] <- counts[chosen]+1L
      if(chosen!=old) changes <- changes+1L
    }
    if(changes==0L) break
  }
  plan
}

# Capped endpoint construction is defined in ocs_cap10_support.R.
