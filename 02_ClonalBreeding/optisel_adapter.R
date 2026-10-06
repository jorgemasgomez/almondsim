.libPaths(c('C:/Users/franc/Documents/Codex/2026-10-04/es/work/r-library',.libPaths()))
optisel_contributions <- function(K,bv,bound=NULL,max_contribution=.5,N=1500) {
 stopifnot(is.matrix(K),nrow(K)==ncol(K),length(bv)==nrow(K),!is.null(rownames(K)),all(is.finite(K)),all(is.finite(bv)))
 phen<-data.frame(Indiv=rownames(K),Sex=NA_character_,BV=bv,isCandidate=TRUE)
 cand<-optiSel::candes(phen,N=N,Kin=K,quiet=TRUE)
 stopifnot(length(max_contribution)%in%c(1L,length(bv)))
 con<-list(ub=setNames(rep_len(max_contribution,length(bv)),rownames(K)))
 method<-if(is.null(bound)) 'min.Kin' else 'max.BV'
 if(!is.null(bound)) con$ub.Kin<-bound
 fit<-optiSel::opticont(method,cand,con,solver='slsqp',quiet=TRUE,maxeval=1000)
 if(!all(fit$info$valid)) stop('optiSel returned contributions violating constraints.')
 fit
}
optisel_compatible_plan <- function(contribution,K,allowed,crosses=20L,max_parents=32L) {
 stopifnot(length(contribution)==nrow(K),all(dim(allowed)==dim(K)))
 contribution <- pmax(0,contribution); contribution<-contribution/sum(contribution)
 desired <- 2*crosses*contribution
 counts <- floor(desired)
 missing <- 2*crosses-sum(counts)
 if(missing>0) {take<-head(order(-(desired-counts),seq_along(counts)),missing);counts[take]<-counts[take]+1L}
 active<-which(counts>0)
 if(length(active)>max_parents) stop('Rounded optiSel contributions require more than 32 parents; no silent truncation.')
 edges<-which(allowed & outer(seq_along(counts)%in%active,seq_along(counts)%in%active,'&'),arr.ind=TRUE)
 if(!nrow(edges)) stop('No compatible mating edges for rounded contributions.')
 A<-sapply(seq_len(nrow(edges)),function(j) as.integer(active==edges[j,1])+as.integer(active==edges[j,2]))
 A<-matrix(A,nrow=length(active))
 fit<-lpSolve::lp('min',objective.in=K[edges],const.mat=A,const.dir=rep('=',length(active)),const.rhs=counts[active],all.int=TRUE,timeout=15L)
 optimal<-fit$status==0L
 if(!optimal) fit<-lpSolve::lp('min',objective.in=rep(0,nrow(edges)),const.mat=A,
   const.dir=rep('=',length(active)),const.rhs=counts[active],all.int=TRUE,timeout=15L)
 if(fit$status!=0) stop('No compatible integer plan for rounded contributions; rerounding or joint optimization is required.')
 plan<-edges[rep(seq_len(nrow(edges)),times=round(fit$solution)),,drop=FALSE]
 realized<-tabulate(as.integer(plan),nbins=length(counts))/(2*crosses)
 stopifnot(nrow(plan)==crosses,all(allowed[plan]),all(tabulate(as.integer(plan),nbins=length(counts))==counts),length(active)<=max_parents)
 list(plan=plan,contributions=realized,continuous_contributions=contribution,mate_allocation_optimal=optimal,
      group_kinship=as.numeric(crossprod(realized,K%*%realized)),
      mate_kinship=mean(K[plan]),rounding_max_error=max(abs(realized-contribution)))
}
# Recovery optimizes INTEGER contributions and compatible pairs together in a
# documented shortlist. It is not optiSel::matings or a global OCS solver.
optisel_recover_plan <- function(contribution,K,bv,allowed,bound,crosses=20L,max_parents=32L,max_cuts=80L) {
 stopifnot(all(is.finite(contribution)),all(is.finite(bv)),is.finite(bound),
           max_parents>=2L,crosses>=1L,isSymmetric(K,tol=1e-8))
 if(min(eigen(K,symmetric=TRUE,only.values=TRUE)$values) < -1e-7)
   stop('Recovery requires a positive semidefinite kinship matrix.')
 allowed<-allowed & !diag(TRUE,nrow(K))
 eligible<-which(rowSums(allowed)+colSums(allowed)>0)
 if(!length(eligible)) stop('Recovery impossible: no compatible pairs.')
 # Keep high-contribution parents and add their best compatible partners.
 ranked<-eligible[order(-contribution[eligible],-bv[eligible],eligible)]
 selected<-integer()
 for(i in ranked) {
   if(i %in% selected) next
   partners<-eligible[allowed[i,eligible]|allowed[eligible,i]]
   partners<-partners[order(-contribution[partners],-bv[partners],partners)]
   add<-if(any(partners %in% selected)) i else unique(c(i,partners[1]))
   if(length(union(selected,add))<=max_parents) selected<-union(selected,add)
   if(length(selected)==max_parents) break
 }
 edges_local<-which(allowed[selected,selected,drop=FALSE],arr.ind=TRUE)
 if(!nrow(edges_local)) stop('Recovery impossible: shortlist has no compatible pairs.')
 edges<-matrix(selected[edges_local],ncol=2)
 B<-matrix(0,length(selected),nrow(edges))
 for(j in seq_len(nrow(edges))) B[,j]<-(as.integer(selected==edges[j,1])+as.integer(selected==edges[j,2]))/(2*crosses)
 Ks<-K[selected,selected,drop=FALSE]
 A<-rbind(rep(1,nrow(edges)),B)
 directions<-c('=',rep('<=',length(selected)))
 rhs<-c(crosses,rep(.5,length(selected)))
 objective<-as.numeric(crossprod(B,bv[selected]))
 # Recovery preserves the optiSel allocation as closely as possible. Continuous
 # deviation variables avoid an unnecessarily difficult gain-only integer search.
 target<-contribution[selected]/sum(contribution[selected])
 nedge<-nrow(edges);nparent<-length(selected)
 A<-cbind(A,matrix(0,nrow(A),nparent))
 A<-rbind(A,cbind(B,-diag(nparent)),cbind(-B,-diag(nparent)))
 directions<-c(directions,rep('<=',2*nparent))
 rhs<-c(rhs,target,-target)
 objective<-c(1e-6*objective,rep(-1,nparent))
 for(iteration in seq_len(max_cuts)) {
   fit<-lpSolve::lp('max',objective,A,directions,rhs,int.vec=seq_len(nedge),timeout=15L)
   counts<-fit$solution[seq_len(nedge)]
   residual<-as.numeric(A%*%fit$solution)
   feasible_incumbent<-all(is.finite(fit$solution)) &&
     max(abs(counts-round(counts)))<1e-7 &&
     abs(sum(counts)-crosses)<1e-7 && all(fit$solution>=-1e-8) &&
     all(residual[directions=='<=']<=rhs[directions=='<=']+1e-7)
   if(!feasible_incumbent) stop('Integer recovery has no usable incumbent (status ',fit$status,'); no invalid plan returned.')
   cnew<-as.numeric(B%*%round(counts))
   kin<-as.numeric(crossprod(cnew,Ks%*%cnew))
   if(kin<=bound+1e-8) {
     plan<-edges[rep(seq_len(nrow(edges)),round(counts)),,drop=FALSE]
     realized<-tabulate(as.integer(plan),nbins=nrow(K))/(2*crosses)
     stopifnot(nrow(plan)==crosses,all(allowed[plan]),sum(realized>0)<=max_parents)
     return(list(plan=plan,contributions=realized,continuous_contributions=contribution,
       group_kinship=kin,mate_kinship=mean(K[plan]),
       rounding_max_error=max(abs(realized-contribution)),recovered=TRUE,
       recovery_method='integer closest-contribution allocation on compatible shortlist with convex kinship cuts',
       recovery_iterations=iteration,recovery_shortlist=rownames(K)[selected],bound=bound,
       recovery_distance_optimal=fit$status==0L))
   }
   # For PSD K, this tangent is a necessary condition for c'Kc <= bound.
   gradient<-2*as.numeric(Ks%*%cnew)
   A<-rbind(A,c(as.numeric(crossprod(B,gradient)),rep(0,nparent)))
   directions<-c(directions,'<=');rhs<-c(rhs,bound+kin)
 }
 stop('Recovery reached its iteration limit; no invalid plan returned.')
}
optisel_validated_plan <- function(K,bv,allowed,bound,crosses=20L,max_parents=32L,max_attempts=12L) {
 continuous_bound<-bound
 for(attempt in seq_len(max_attempts)) {
   fit<-optisel_contributions(K,bv,continuous_bound)
   c<-fit$parent$oc[match(rownames(K),fit$parent$Indiv)]
   plan<-tryCatch(optisel_compatible_plan(c,K,allowed,crosses,max_parents),error=function(e)e)
   if(inherits(plan,'error')) {
     recovered<-optisel_recover_plan(c,K,bv,allowed,bound,crosses,max_parents)
     recovered$fit<-fit;recovered$recovery_reason<-conditionMessage(plan)
     recovered$rounding_attempts<-attempt
     return(recovered)
   }
   if(plan$group_kinship<=bound+1e-8) {
     plan$fit<-fit;plan$bound<-bound;plan$continuous_bound<-continuous_bound
     plan$rounding_attempts<-attempt
     plan$recovered<-FALSE
     return(plan)
   }
   continuous_bound <- continuous_bound-max(2*(plan$group_kinship-bound),1e-4)
 }
 recovered<-optisel_recover_plan(c,K,bv,allowed,bound,crosses,max_parents)
 recovered$fit<-fit;recovered$recovery_reason<-'Rounding exceeded the kinship bound after adjustment.'
 recovered$rounding_attempts<-max_attempts
 recovered
}
