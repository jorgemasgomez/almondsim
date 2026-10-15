# Native optiSel contribution -> offspring count -> mating workflow.
# This adapter defines almond biology and validates results; it does not repair,
# locally refine or substitute contribution vectors returned by optiSel.
if(!requireNamespace('optiSel',quietly=TRUE)) {
 lib<-Sys.getenv('ALMONDSIM_OCS_R_LIBRARY','C:/Users/franc/Documents/Codex/2026-10-04/es/work/r-library')
 if(dir.exists(lib)).libPaths(c(lib,.libPaths()))
}
if(!requireNamespace('optiSel',quietly=TRUE)||!requireNamespace('lpSolve',quietly=TRUE))
 stop('OCS requires optiSel and lpSolve.')

optisel_offspring_context <- function(K,bv,allowed,N,max_contribution=.10) {
 stopifnot(is.matrix(K),nrow(K)==ncol(K),length(bv)==nrow(K),
   identical(dim(allowed),dim(K)),!anyNA(allowed),all(is.finite(K)),all(is.finite(bv)),
   length(N)==1L,is.finite(N),N>=1,N==floor(N),
   max_contribution>0,max_contribution<=1,max(abs(K-t(K)))<1e-8)
 diag(allowed)<-FALSE
 mothers<-which(rowSums(allowed)>0);fathers<-which(colSums(allowed)>0)
 if(!length(mothers)||!length(fathers))stop('No eligible maternal or paternal candidates.')
 index<-c(mothers,fathers)
 ids<-c(paste0('M_',mothers),paste0('P_',fathers))
 Kr<-K[index,index,drop=FALSE];dimnames(Kr)<-list(ids,ids)
 scaled<-if(sd(bv)>0) as.numeric(scale(bv)) else rep(0,length(bv))
 phen<-data.frame(Indiv=ids,Sex=c(rep('female',length(mothers)),rep('male',length(fathers))),
   BV=scaled[index],isCandidate=TRUE)
 # Native linear trait constraints link both role copies of the same tree.
 # Thus ub=0.10 applies to the combined tree, not 0.10 separately per role.
 con<-list(ub=setNames(rep(max_contribution,length(index)),ids))
 for(i in intersect(mothers,fathers)) {
   trait<-paste0('TreeCap',i)
   phen[[trait]]<-as.numeric(index==i)
   con[[paste0('ub.',trait)]]<-max_contribution
 }
 cand<-optiSel::candes(phen,N=N,Kin=Kr,quiet=TRUE)
 q<-cand$kinship$Kin
 if(min(eigen(K,symmetric=TRUE,only.values=TRUE)$values)< -1e-8)
   stop('OCS requires a positive-semidefinite genomic coancestry matrix.')
 if(max(abs(q$a))>1e-10)stop('Unexpected linear term in optiSel coancestry representation.')
 list(K=K,bv=bv,allowed=allowed,N=as.integer(N),cap=max_contribution,
   cand=cand,con=con,drift=q$d,index=index,ids=ids)
}

optisel_native_fit <- function(context,method,bound=NULL) {
 con<-context$con
 if(!is.null(bound))con$ub.Kin<-bound+context$drift
 fit<-optiSel::opticont(method,context$cand,con,solver='slsqp',quiet=TRUE,maxeval=3000)
 if(!all(fit$info$valid))stop('Native optiSel continuous constraints failed.')
 oc<-fit$parent$oc[match(context$ids,fit$parent$Indiv)]
 if(anyNA(oc)||any(!is.finite(oc))||any(oc< -1e-7))stop('Invalid optiSel contribution vector.')
 combined<-numeric(nrow(context$K))
 for(j in seq_along(oc))combined[context$index[j]]<-combined[context$index[j]]+oc[j]
 if(abs(sum(combined)-1)>1e-6||max(combined)>context$cap+1e-7)
   stop('Native optiSel combined individual contribution cap failed.')
 kin<-as.numeric(crossprod(combined,context$K%*%combined))
 if(!is.null(bound)&&kin>bound+1e-7)stop('Native continuous coancestry bound failed.')
 list(fit=fit,contribution=combined,kinship=kin,gain=sum(combined*context$bv))
}

optisel_matings_offspring <- function(phen,K,allowed,index_by_id) {
 phen<-phen[phen$n>0,,drop=FALSE]
 # Native matings() retains a matrix when at least two candidates have each role.
 # A valid combined 10% cap necessarily needs at least five of each role.
 if(sum(phen$Sex=='female')<2L||sum(phen$Sex=='male')<2L)
   stop('Too few contributing parental roles for native matings().')
 checked<-getFromNamespace('checkphen','optiSel')(phen,columns=c('Indiv','Sex'),quiet=TRUE,na.Sex=FALSE)
 sire<-checked$Indiv[checked$Sex=='male'];dam<-checked$Indiv[checked$Sex=='female']
 cost<-K[index_by_id[sire],index_by_id[dam],drop=FALSE]
 dimnames(cost)<-list(sire,dam)
 valid<-t(allowed[index_by_id[dam],index_by_id[sire],drop=FALSE])
 forbidden<-which(!as.vector(valid));solver_status<-NA_integer_
 # Use matings()'s public solver interface solely to impose forbidden edges.
 # Its native objective and maternal/paternal offspring margins are preserved.
 solver_hook<-function(obj,mat,dir,rhs,types,bounds,max,...) {
   stopifnot(length(obj)==length(valid),max(abs(obj-as.vector(cost)))<1e-12,is.null(bounds))
   if(length(forbidden)) {
     extra<-matrix(0,length(forbidden),length(obj))
     extra[cbind(seq_along(forbidden),forbidden)]<-1
     mat<-rbind(mat,extra);dir<-c(dir,rep('==',length(forbidden)));rhs<-c(rhs,rep(0,length(forbidden)))
   }
   z<-lpSolve::lp(if(max)'max'else'min',obj,mat,ifelse(dir=='==','=',dir),rhs,
     all.int=TRUE,timeout=30L)
   solver_status<<-z$status
   if(z$status!=0L)stop('Native matings could not allocate the fixed offspring counts; lpSolve status ',z$status)
   stopifnot(max(abs(z$solution-round(z$solution)))<1e-7)
   list(solution=round(z$solution),status=0L)
 }
 mating<-optiSel::matings(phen,cost,solver=solver_hook)
 plan<-cbind(unname(index_by_id[mating$Dam]),unname(index_by_id[mating$Sire]))
 weights<-as.integer(mating$n)
 if(anyNA(plan)||any(weights<1L)||!all(allowed[plan]))stop('Invalid native mating allocation.')
 wanted<-setNames(phen$n,phen$Indiv)
 for(sex in c('female','male')) {
   names_role<-if(sex=='female')mating$Dam else mating$Sire
   got<-tapply(weights,names_role,sum)
   expected<-wanted[phen$Indiv[phen$Sex==sex]]
   if(anyNA(got[names(expected)])||any(got[names(expected)]!=expected))
     stop('matings() changed the offspring counts supplied by noffspring().')
 }
 list(plan=plan,family_sizes=weights,matings=mating,
   solver_status=solver_status,forbidden_edges=length(forbidden))
}

optisel_offspring_profile <- function(K,bv,allowed,N,alpha,max_contribution=.10) {
 stopifnot(length(alpha)==1L,is.finite(alpha),alpha>=0,alpha<=1)
 context<-optisel_offspring_context(K,bv,allowed,N,max_contribution)
 diversity<-optisel_native_fit(context,'min.Kin')
 gain<-optisel_native_fit(context,'max.BV')
 if(gain$kinship<diversity$kinship-1e-7)stop('Inconsistent native continuous coancestry endpoints.')
 low<-diversity$kinship;high<-max(low,gain$kinship)
 B<-low+alpha*(high-low)
 optimum<-optisel_native_fit(context,'max.BV',B)
 phen<-optimum$fit$parent
 counts<-optiSel::noffspring(phen,N=N,random=FALSE)
 phen$n<-counts$nOff[match(phen$Indiv,counts$Indiv)]
 if(anyNA(phen$n)||any(phen$n<0)||any(phen$n!=floor(phen$n))||
   sum(phen$n[phen$Sex=='female'])!=N||sum(phen$n[phen$Sex=='male'])!=N)
   stop('Native noffspring() did not preserve the offspring budget.')
 index_by_id<-setNames(context$index,context$ids)
 rounded<-numeric(nrow(K))
 for(j in seq_len(nrow(phen)))rounded[index_by_id[phen$Indiv[j]]]<-
   rounded[index_by_id[phen$Indiv[j]]]+phen$n[j]/(2*N)
 rounded_kin<-as.numeric(crossprod(rounded,K%*%rounded))
 # Preserve native counts. Individual cap/biology remain hard constraints;
 # coancestry departures introduced by rounding are reported without stopping.
 if(max(rounded)>max_contribution+1e-10)
   stop('Native noffspring rounding exceeded the combined 10% cap; no custom repair was applied.')
 kinship_excess<-max(0,rounded_kin-B)
 kinship_excess_relative<-if(B>0)kinship_excess/B else if(kinship_excess==0)0 else Inf
 if(kinship_excess>1e-8)
   warning(sprintf('Native noffspring rounding: coancestry %.12g exceeds target %.12g by %.12g (%.6g%%); preserving native offspring counts.',
     rounded_kin,B,kinship_excess,100*kinship_excess_relative),call.=FALSE)
 final<-optisel_matings_offspring(phen,K,context$allowed,index_by_id)
 met<-ocs_metrics(final$plan,K,bv,weights=final$family_sizes)
 stopifnot(sum(final$family_sizes)==N,max(met$contribution)<=max_contribution+1e-10,
   max(abs(met$contribution-rounded))<1e-10)
 list(plan=final$plan,family_sizes=final$family_sizes,metrics=met,bound=B,minK=low,maxK=high,alpha=alpha,
   audit=list(implementation='optiSel::opticont -> optiSel::noffspring -> optiSel::matings',
     continuous_gain=optimum$gain,continuous_kinship=optimum$kinship,
     continuous_contributions=optimum$contribution,continuous_fit=optimum$fit,rounded_contributions=rounded,
     offspring_counts=phen[c('Indiv','Sex','n')],role_to_tree=index_by_id,
     endpoint_diversity=diversity,endpoint_gain=gain,
     final_matings=final$matings,forbidden_edges=final$forbidden_edges,
     solver_status=final$solver_status,max_individual_contribution=max_contribution,
     offspring=N,family_sizes=final$family_sizes,distinct_pairs=nrow(final$plan),
     continuous_bound=B,custom_rounding_allowance=FALSE,
     rounded_bound_satisfied=kinship_excess<=1e-8,
     rounded_kinship_excess=kinship_excess,rounded_kinship_excess_relative=kinship_excess_relative,
     max_parents=NA_integer_,custom_contribution_refinement=FALSE))
}
