# Allocation tests only: no breeding generations, founder bank or ASReml needed.
# Run from the repository root: Rscript 02_ClonalBreeding/tests/check_ocs_cap10.R
root<-normalizePath('02_ClonalBreeding',winslash='/',mustWork=TRUE)
options(almond.input_dir=root)
source(file.path(root,'redesigned_schemes.R'))
source(file.path(root,'ocs_discrete.R'))
source(file.path(root,'manual_ocs_adapter.R'))
stopifnot(redesign_config()$design_version=='PT_AT_ET_optiSel_cap10_v8',redesign_config()$ocs_max_contribution==.10)
K<-diag(3);allowed<-matrix(TRUE,3,3);diag(allowed)<-FALSE
for(plan in list(cbind(rep(1,20),rep(2,20)),cbind(rep(c(1,2),each=10),rep(3,20)),cbind(rep(1,20),rep(c(2,3),each=10)))) {
 z<-optisel_matings_exact(plan,K,allowed)
 stopifnot(nrow(z$plan)==20L,z$solver_status==0L)
}
K<-diag(24);dimnames(K)<-list(as.character(1:24),as.character(1:24));bv<-seq_len(24)
for(kind in c('normal','speed')) {
 allowed<-matrix(TRUE,24,24);diag(allowed)<-FALSE
 if(kind=='speed')allowed[9:24,]<-FALSE
 for(a in c(.10,.20,.40)) {
  z<-optisel_matings_profile(K,bv,allowed,20,32,a,seed=1901L)
  stopifnot(max(tabulate(z$plan,nbins=24))<=4L,length(unique(as.integer(z$plan)))>=10L,
   length(unique(z$plan[,1]))>=5L,length(unique(z$plan[,2]))>=5L,
   all(allowed[z$plan]),z$metrics$kinship<=z$bound+1e-8,z$audit$solver_status==0L,
   max(z$audit$continuous_contributions)<=.10+1e-7)
 }
}
bad<-matrix(FALSE,24,24);bad[1:4,5:24]<-TRUE
stopifnot(inherits(tryCatch(ocs_cap_transport(K,bv,bad,20,32),error=function(e)e),'error'))
cat('PASS: singleton dimensions, continuous/integer 10% cap, both roles, all profiles, impossible maternal budget.\n')
