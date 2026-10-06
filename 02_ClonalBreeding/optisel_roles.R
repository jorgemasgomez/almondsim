# Standard sex-specific optiSel contribution optimization. The same tree may
# occur in both roles; role IDs are distinct but its genotype/kinship is shared.
optisel_role_fit <- function(K,bv,maternal,paternal,bound=NULL,upper=NULL) {
 m<-which(maternal);p<-which(paternal);index<-c(m,p)
 ids<-c(paste0(rownames(K)[m],'_M'),paste0(rownames(K)[p],'_P'))
 Kr<-K[index,index,drop=FALSE];dimnames(Kr)<-list(ids,ids)
 phen<-data.frame(Indiv=ids,Sex=c(rep('female',length(m)),rep('male',length(p))),
   BV=bv[index],isCandidate=TRUE)
 candidates<-optiSel::candes(phen,N=1500,Kin=Kr,quiet=TRUE)
 if(is.null(upper)) upper<-rep(.5,length(index))
 if(all(maternal==paternal)) {
   selected<-m
   cap<-pmin(.5,2*pmin(upper[seq_along(m)],upper[length(m)+seq_along(p)]))
   fit<-optisel_contributions(K[selected,selected,drop=FALSE],bv[selected],bound,cap)
   c<-fit$parent$oc[match(rownames(K)[selected],fit$parent$Indiv)]
   return(list(fit=fit,contributions=c(c/2,c/2),index=index,mothers=seq_along(m),
     fathers=length(m)+seq_along(p),upper=upper))
 }
 constraints<-list(ub=setNames(upper,ids))
 if(!is.null(bound)) constraints$ub.Kin<-bound
 fit<-optiSel::opticont(if(is.null(bound)) 'min.Kin' else 'max.BV',candidates,
   constraints,solver='slsqp',quiet=TRUE,maxeval=1500)
 if(!all(fit$info$valid)) stop('Role-specific optiSel constraints were not satisfied.')
 c<-fit$parent$oc[match(ids,fit$parent$Indiv)]
 list(fit=fit,contributions=pmax(0,c),index=index,mothers=seq_along(m),
      fathers=length(m)+seq_along(p),upper=upper)
}
optisel_transport <- function(role,K,allowed,crosses=20L) {
 round_role<-function(indices) {
   desired<-crosses*role$contributions[indices]/sum(role$contributions[indices])
   n<-floor(desired);missing<-crosses-sum(n)
   if(missing>0) {chosen<-head(order(-(desired-n),seq_along(n)),missing);n[chosen]<-n[chosen]+1L}
   n
 }
 mother_counts<-round_role(role$mothers);father_counts<-round_role(role$fathers)
 mothers<-role$index[role$mothers];fathers<-role$index[role$fathers]
 if(identical(mothers,fathers)) {
   desired<-2*crosses*(role$contributions[role$mothers]+role$contributions[role$fathers])
   counts<-floor(desired);missing<-2*crosses-sum(counts)
   if(missing>0) {chosen<-head(order(-(desired-counts),seq_along(counts)),missing);counts[chosen]<-counts[chosen]+1L}
   mother_counts<-floor(counts/2)
   odd<-which(counts%%2==1)
   needed<-crosses-sum(mother_counts)
   if(needed>0) {
     order_odd<-odd[order(-(2*crosses*role$contributions[role$mothers][odd]-mother_counts[odd]),odd)]
     mother_counts[head(order_odd,needed)]<-mother_counts[head(order_odd,needed)]+1L
   }
   father_counts<-counts-mother_counts
 }
 active_m<-which(mother_counts>0);active_p<-which(father_counts>0)
 valid<-allowed[mothers[active_m],fathers[active_p],drop=FALSE]
 edges<-which(valid,arr.ind=TRUE)
 if(!nrow(edges)) stop('No compatible edges for role-specific contributions.')
 A<-matrix(0,length(active_m)+length(active_p),nrow(edges))
 for(j in seq_len(nrow(edges))) {A[edges[j,1],j]<-1;A[length(active_m)+edges[j,2],j]<-1}
 rhs<-c(mother_counts[active_m],father_counts[active_p])
 pair<-cbind(mothers[active_m][edges[,1]],fathers[active_p][edges[,2]])
 # A bipartite transport incidence matrix is totally unimodular: integer
 # supplies/demands yield integer vertex solutions without branch-and-bound.
 fit<-lpSolve::lp('min',K[pair],A,rep('=',length(rhs)),rhs,timeout=15L)
 if(fit$status!=0L) stop('Role-specific compatible transport has no certified solution.')
 stopifnot(max(abs(fit$solution-round(fit$solution)))<1e-7)
 plan<-pair[rep(seq_len(nrow(pair)),round(fit$solution)),,drop=FALSE]
 c<-tabulate(as.integer(plan),nbins=nrow(K))/(2*crosses)
 stopifnot(nrow(plan)==crosses,all(allowed[plan]),
   all(tabulate(plan[,1],nbins=nrow(K))[mothers]==mother_counts),
   all(tabulate(plan[,2],nbins=nrow(K))[fathers]==father_counts))
 list(plan=plan,contributions=c,group_kinship=as.numeric(crossprod(c,K%*%c)),
   mother_counts=mother_counts,father_counts=father_counts,mate_allocation_optimal=TRUE)
}
optisel_role_minimum <- function(K,bv,allowed,maternal,paternal,crosses,max_parents) {
 upper<-rep(.5,sum(maternal)+sum(paternal));shortlist<-seq_len(nrow(K))
 for(attempt in 1:30) {
   role<-optisel_role_fit(K,bv,maternal,paternal,upper=upper)
   plan<-tryCatch(optisel_transport(role,K,allowed,crosses),error=function(e)e)
   if(inherits(plan,'error')) {
     # Keep enough role capacity to fill each half of the mating budget.
     changed<-FALSE
     for(indices in list(role$fathers,role$mothers)) {
       order<-indices[order(-role$contributions[indices],indices)]
       for(j in order) if(upper[j]>.025 && sum(upper[indices])-.025>=.5-1e-8) {
         upper[j]<-upper[j]-.025;changed<-TRUE;break
       }
     }
     if(!changed) stop('No compatible role-specific diversity endpoint could be constructed.')
     next
   }
   if(sum(plan$contributions>0)>max_parents) {
     combined<-numeric(nrow(K))
     totals<-rowsum(role$contributions,role$index,reorder=FALSE)
     combined[as.integer(rownames(totals))]<-totals[,1]
     shortlist<-order(-combined,seq_along(combined))[seq_len(max_parents)]
     maternal<-maternal & seq_along(maternal)%in%shortlist
     paternal<-paternal & seq_along(paternal)%in%shortlist
     upper<-rep(.5,sum(maternal)+sum(paternal));next
   }
   plan$role<-role;plan$shortlist<-shortlist
   return(plan)
 }
 stop('Could not construct a feasible diversity endpoint within 30 attempts.')
}
optisel_role_solution <- function(K,bv,allowed,maternal,paternal,alpha,crosses=20L,max_parents=32L) {
 minimum<-optisel_role_minimum(K,bv,allowed,maternal,paternal,crosses,max_parents)
 pairs<-which(allowed,arr.ind=TRUE)
 if(!nrow(pairs)) stop('No eligible S-compatible parental pair.')
 best<-pairs[which.max((bv[pairs[,1]]+bv[pairs[,2]])/2),]
 endpoint<-tabulate(best,nbins=nrow(K))/2
 maxK<-as.numeric(crossprod(endpoint,K%*%endpoint))
 minK<-minimum$group_kinship
 collapsed<-minK>maxK
 if(collapsed) {
   minimum$plan<-matrix(rep(best,each=crosses),ncol=2)
   minimum$contributions<-endpoint;minimum$group_kinship<-maxK
   minK<-maxK
 }
 bound<-minK+alpha*(maxK-minK)
 upper<-rep(.5,sum(maternal)+sum(paternal));original_bound<-bound
 reasons<-character();shortlist<-seq_len(nrow(K));fit_bound<-bound
 for(attempt in 1:30) {
   role<-tryCatch(optisel_role_fit(K,bv,maternal,paternal,fit_bound,upper),error=function(e)e)
   if(inherits(role,'error')) {reasons<-c(reasons,conditionMessage(role));break}
   ccombined<-rowsum(role$contributions,group=role$index,reorder=FALSE)
   combined<-numeric(nrow(K));combined[as.integer(rownames(ccombined))]<-ccombined[,1]
   # Selfing is prohibited; no tree can supply more than half the total.
   dominant<-which(combined>.5+1e-7)
   if(length(dominant)) {
     for(i in dominant) {
       j<-intersect(which(role$index==i),role$fathers)
       other<-sum(role$contributions[role$index==i])-sum(role$contributions[j])
       if(length(j)) upper[j]<-pmin(upper[j],pmax(0,.5-other))
     }
     reasons<-c(reasons,'Individual contribution above 50%; role caps reduced.')
     next
   }
   plan<-tryCatch(optisel_transport(role,K,allowed,crosses),error=function(e)e)
   if(inherits(plan,'error')) {
     # Restrict highly concentrated role contributions; retain both role budgets.
     for(indices in list(role$fathers,role$mothers)) {
       order<-indices[order(-role$contributions[indices],indices)]
       for(j in order) if(upper[j]>.025 && sum(upper[indices])-.025>=.5-1e-8) {
         upper[j]<-upper[j]-.025;break
       }
     }
     reasons<-c(reasons,conditionMessage(plan));next
   }
   if(sum(plan$contributions>0)>max_parents) {
     shortlist<-order(-combined,seq_along(combined))[seq_len(max_parents)]
     if(!any(maternal[shortlist])||!any(paternal[shortlist])) stop('Shortlist lost a required parental role.')
     maternal<-maternal & seq_along(maternal)%in%shortlist
     paternal<-paternal & seq_along(paternal)%in%shortlist
     upper<-rep(.5,sum(maternal)+sum(paternal))
     reasons<-c(reasons,'Reoptimized contributions within a 32-parent shortlist.');next
   }
   if(plan$group_kinship<=original_bound+1e-8) {
     plan$fit<-role$fit;plan$role<-role
     plan$bound<-original_bound;plan$minK<-minK;plan$maxK<-maxK;plan$alpha<-alpha
     plan$attempts<-attempt;plan$recovery_reasons<-reasons
     plan$fallback_to_diversity_endpoint<-FALSE;plan$endpoints_collapsed<-collapsed
     plan$shortlist<-rownames(K)[shortlist]
     return(plan)
   }
   fit_bound<-fit_bound-max(2*(plan$group_kinship-original_bound),1e-4)
   reasons<-c(reasons,'Kinship bound tightened after role-specific rounding.')
 }
 # A feasible endpoint is retained before gain optimization. A solver failure
 # can use it without relaxing the realized bound or creating illegal crosses.
 minimum$bound<-original_bound;minimum$minK<-minK;minimum$maxK<-maxK;minimum$alpha<-alpha
 minimum$attempts<-attempt;minimum$recovery_reasons<-reasons
 minimum$fallback_to_diversity_endpoint<-TRUE;minimum$endpoints_collapsed<-collapsed
 stopifnot(minimum$group_kinship<=original_bound+1e-8,all(allowed[minimum$plan]),
   sum(minimum$contributions>0)<=max_parents,nrow(minimum$plan)==crosses)
 minimum
}
