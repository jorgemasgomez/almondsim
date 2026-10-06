# Genomic OCS on realized integer mating plans. Heuristic local search;
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

ocs_endpoints <- function(K,bv,allowed,crosses,max_parents,restarts=4L,sweeps=12L) {
  pairs <- which(allowed,arr.ind=TRUE)
  if(!nrow(pairs)) stop('No S-compatible non-self crosses in the OCS candidate pool.')
  pair_gain <- (bv[pairs[,1]]+bv[pairs[,2]])/2
  best <- which.max(pair_gain)
  gain_plan <- matrix(rep(pairs[best,],each=crosses),nrow=crosses,ncol=2)
  # Repeating a family is permitted in both comparison arms; budget counts
  # crossing events rather than requiring 20 distinct parental pairs.
  minimum <- gain_plan
  bestK <- ocs_metrics(minimum,K,bv)$kinship
  for(restart in seq_len(restarts)) {
    pool <- sample.int(nrow(K),min(max_parents,nrow(K)))
    feasible <- pairs[pairs[,1]%in%pool & pairs[,2]%in%pool,,drop=FALSE]
    if(!nrow(feasible)) next
    plan <- feasible[sample.int(nrow(feasible),crosses,replace=TRUE),,drop=FALSE]
    plan <- ocs_improve(plan,K,bv,allowed,max_parents,'diversity',sweeps=sweeps)
    score <- ocs_metrics(plan,K,bv)$kinship
    if(score<bestK) {minimum <- plan;bestK <- score}
  }
  list(diversity=minimum,gain=gain_plan,
       minK=bestK,maxK=ocs_metrics(gain_plan,K,bv)$kinship)
}

ocs_profile_plan <- function(K,bv,allowed,crosses,max_parents,alpha,
                             restarts=4L,sweeps=12L,seed=1L) {
  stopifnot(alpha>=0,alpha<=1,max_parents>=2,crosses>0)
  # Optimization randomness must not consume the breeding simulation stream.
  had_seed <- exists('.Random.seed',.GlobalEnv,inherits=FALSE)
  if(had_seed) saved <- get('.Random.seed',.GlobalEnv)
  on.exit(if(had_seed) assign('.Random.seed',saved,.GlobalEnv) else
            if(exists('.Random.seed',.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv))
  set.seed(seed)
  endpoints <- ocs_endpoints(K,bv,allowed,crosses,max_parents,restarts,sweeps)
  bound <- endpoints$minK+alpha*(endpoints$maxK-endpoints$minK)
  plan <- endpoints$diversity
  steps <- sort(unique(c(.05,.25,.6,alpha)))
  for(step in steps[steps<=alpha]) {
    step_bound <- endpoints$minK+step*(endpoints$maxK-endpoints$minK)
    plan <- ocs_improve(plan,K,bv,allowed,max_parents,'gain',step_bound,sweeps)
  }
  # Optimize mate allocation while keeping the genetic contributions unchanged.
  # Swapping fathers preserves every count and therefore the OCS bound exactly.
  for(pass in 1:3) for(i in seq_len(crosses)) for(j in seq_len(crosses)) {
    if(i>=j) next
    a <- plan[i,1];b <- plan[i,2];c <- plan[j,1];d <- plan[j,2]
    if(allowed[a,d] && allowed[c,b] && K[a,d]+K[c,b]<K[a,b]+K[c,d]-1e-12) {
      plan[i,2] <- d;plan[j,2] <- b
    }
  }
  metrics <- ocs_metrics(plan,K,bv)
  stopifnot(nrow(plan)==crosses,metrics$parents<=max_parents,
            all(allowed[plan]),metrics$kinship<=bound+1e-8)
  list(plan=plan,metrics=metrics,bound=bound,minK=endpoints$minK,maxK=endpoints$maxK,alpha=alpha)
}

redesign_ocs <- function(state,candidates,model) {
  if(is.null(state$ocs_frequency)) stop('Missing fixed founder SNP frequencies for OCS.')
  pool <- redesign_merge_latest(list(state$parents,candidates))
  pool <- redesign_predict(pool,'GS',model)
  K <- ocs_kinship(pool,state$ocs_frequency)
  allowed <- ocs_compatibility(pool,locuscompt_position_vector,self_compatible_allele)
  solution <- ocs_profile_plan(K,as.numeric(ebv(pool)),allowed,state$config$crosses,
                              state$config$parents,state$config$ocs_alpha,
                              restarts=state$config$ocs_restarts,sweeps=state$config$ocs_sweeps,
                              seed=state$config$seed+as.integer(state$year))
  active <- which(solution$metrics$counts>0)
  selected <- pool[active]
  old_ids <- state$parents@id
  solution$new <- sum(!selected@id %in% old_ids)
  solution$weighted_age <- sum(solution$metrics$contribution*
                               (state$year-state$birth[as.character(pool@id)]))
  solution$contributions <- data.frame(id=as.character(pool@id[active]),
                                       contribution=solution$metrics$contribution[active])
  solution$cross_ids <- matrix(as.character(pool@id[solution$plan]),ncol=2)
  solution$parents <- selected
  solution
}
