# Loaded by the OCS backend; keeps the installed optiSel package unchanged.
ocs_cap_slots<-function(crosses) floor(2*crosses*.10+1e-8)
optisel_matings_dimension_safe<-local({
 f<-optiSel::matings
 code<-deparse(f)
 stopifnot(sum(grepl('Kin <- Kin[Sire, Dam]',code,fixed=TRUE))==1L)
 code<-sub('Kin <- Kin[Sire, Dam]','Kin <- Kin[Sire, Dam, drop = FALSE]',code,fixed=TRUE)
 f<-eval(parse(text=code),envir=asNamespace('optiSel'))
 f
})

# Feasible capped integer construction on a bounded candidate shortlist.
# Each individual has one combined maternal+paternal participation budget.
ocs_cap_transport<-function(K,bv,allowed,crosses,max_parents,mode='gain',role=NULL,timeout=15L) {
 n<-length(bv);cap<-ocs_cap_slots(crosses)
 m<-which(rowSums(allowed)>0);p<-which(colSums(allowed)>0)
 if(length(m)*cap<crosses||length(p)*cap<crosses||n*cap<2*crosses||max_parents*cap<2*crosses)
   stop('Insufficient eligible parents under the 10% cap; biological constraints cannot be relaxed.')
 rank_m<-m[order(-bv[m])];rank_p<-p[order(-bv[p])]
 if(is.null(role)) {
   pool<-unique(c(head(rank_m,24),head(rank_p,24),head(order(diag(K)),16)))
 } else {
   combined<-numeric(n);t<-rowsum(role$contributions,role$index,reorder=FALSE);combined[as.integer(rownames(t))]<-t[,1]
   pool<-unique(c(head(order(-combined),48),head(rank_m,16),head(rank_p,16)))
 }
 # Expand shortlist if necessary rather than return an uncapped backup.
 for(extra in c(0L,24L,48L)) {
   keep<-unique(c(pool,head(rank_m,24+extra),head(rank_p,24+extra)))
   edges<-which(allowed[keep,keep,drop=FALSE],arr.ind=TRUE)
   if(!nrow(edges))next
   P<-length(keep);E<-nrow(edges);pairs<-matrix(keep[edges],ncol=2)
   inc<-matrix(0,P,E)
   for(j in seq_len(E)) {inc[edges[j,1],j]<-1;inc[edges[j,2],j]<-inc[edges[j,2],j]+1}
   A<-matrix(0,2L+P,E+P)
   A[1,seq_len(E)]<-1;A[2,E+seq_len(P)]<-1
   A[2+seq_len(P),seq_len(E)]<-inc
   A[cbind(2+seq_len(P),E+seq_len(P))]<- -cap
   rhs<-c(crosses,max_parents,rep(0,P));dirs<-c('=','<=',rep('<=',P))
   cost<-if(mode=='gain') -(bv[pairs[,1]]+bv[pairs[,2]])/2 else K[pairs]
   if(!is.null(role)) cost<-cost - (combined[pairs[,1]]+combined[pairs[,2]])*max(1,sd(cost))
   cost<-cost/max(1,max(abs(cost)))
   fit<-lpSolve::lp('min',c(cost,rep(0,P)),A,dirs,rhs,int.vec=seq_len(E),binary.vec=E+seq_len(P),timeout=timeout)
   if(fit$status!=0L)next
   counts<-as.integer(round(fit$solution[seq_len(E)]))
   plan<-pairs[rep(seq_len(E),counts),,drop=FALSE]
   stopifnot(nrow(plan)==crosses,all(allowed[plan]),max(tabulate(plan,nbins=n))<=cap,
      length(unique(as.integer(plan)))<=max_parents,all(plan[,1]!=plan[,2]))
   return(plan)
 }
 stop('No certified feasible capped integer allocation found on the recovery shortlist.')
}

ocs_endpoints<-function(K,bv,allowed,crosses,max_parents,restarts=4L,sweeps=12L) {
 gain_plan<-ocs_cap_transport(K,bv,allowed,crosses,max_parents,'gain')
 gain_plan<-ocs_improve(gain_plan,K,bv,allowed,max_parents,'gain',Inf,sweeps)
 minimum<-gain_plan;bestK<-ocs_metrics(minimum,K,bv)$kinship
 for(i in seq_len(restarts)) {
  seed<-if(i==1L) gain_plan else ocs_cap_transport(K,bv,allowed,crosses,max_parents,'diversity')
  plan<-ocs_improve(seed,K,bv,allowed,max_parents,'diversity',sweeps=sweeps)
  score<-ocs_metrics(plan,K,bv)$kinship
  if(score<bestK){minimum<-plan;bestK<-score}
 }
 list(diversity=minimum,gain=gain_plan,minK=bestK,maxK=ocs_metrics(gain_plan,K,bv)$kinship)
}

optisel_recover_transport<-function(role,K,allowed,crosses,max_parents,bv,bound=NULL,timeout=2L) {
 plan<-ocs_cap_transport(K,bv,allowed,crosses,max_parents,'gain',role=role,timeout=max(timeout,15L))
 list(plan=plan,recovery_stage='capped_integer_role_recovery')
}
