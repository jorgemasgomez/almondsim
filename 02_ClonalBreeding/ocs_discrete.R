# Genomic coancestry, biological compatibility and offspring-weighted metrics.

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

# Weights are offspring numbers for variable-size families, or one per event.
ocs_metrics <- function(plan,K,bv,weights=rep(1,nrow(plan))) {
 stopifnot(ncol(plan)==2L,length(weights)==nrow(plan),!anyNA(plan),
   all(is.finite(weights)),all(weights>0))
 counts<-numeric(nrow(K))
 for(j in seq_len(nrow(plan)))counts[plan[j,]]<-counts[plan[j,]]+weights[j]
 contribution<-counts/(2*sum(weights));active<-which(counts>0)
 list(kinship=as.numeric(crossprod(contribution[active],K[active,active,drop=FALSE]%*%contribution[active])),
   gain=sum(contribution*bv),parents=length(active),counts=counts,contribution=contribution,
   mate_kinship=sum(weights*K[plan])/sum(weights))
}
