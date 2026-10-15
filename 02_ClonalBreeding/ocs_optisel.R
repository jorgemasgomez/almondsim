if(!requireNamespace('optiSel',quietly=TRUE)) {
  local_library<-Sys.getenv('ALMONDSIM_OCS_R_LIBRARY','C:/Users/franc/Documents/Codex/2026-10-04/es/work/r-library')
  if(dir.exists(local_library)).libPaths(c(local_library,.libPaths()))
}
if(!requireNamespace('optiSel',quietly=TRUE)||!requireNamespace('lpSolve',quietly=TRUE))stop('OCS requires optiSel and lpSolve; configure ALMONDSIM_OCS_R_LIBRARY.')
# optiSel 2.1.0: exact forbidden-edge constraints through its solver hook.
optisel_matings_exact <- function(plan,K,allowed) {
  n<-nrow(K);crosses<-nrow(plan)
  mothers<-which(tabulate(plan[,1],n)>0);fathers<-which(tabulate(plan[,2],n)>0)
  idsM<-paste0('M_',mothers);idsP<-paste0('P_',fathers)
  phen<-data.frame(Indiv=c(idsM,idsP),Sex=c(rep('female',length(mothers)),rep('male',length(fathers))),
    n=c(tabulate(plan[,1],n)[mothers],tabulate(plan[,2],n)[fathers]))
  # Match the exact preprocessing and column-major variable order in matings().
  checked<-getFromNamespace('checkphen','optiSel')(phen,columns=c('Indiv','Sex'),quiet=TRUE,na.Sex=FALSE)
  sire<-checked$Indiv[checked$Sex=='male'];dam<-checked$Indiv[checked$Sex=='female']
  index<-setNames(c(mothers,fathers),c(idsM,idsP))
  cost<-K[index[sire],index[dam],drop=FALSE];dimnames(cost)<-list(sire,dam)
  valid<-t(allowed[index[dam],index[sire],drop=FALSE])
  forbidden<-which(!as.vector(valid));solver_status<-NA_integer_
  solver_hook<-function(obj,mat,dir,rhs,types,bounds,max,...) {
    stopifnot(length(obj)==length(valid),max(abs(obj-as.vector(cost)))<1e-12)
    if(length(forbidden)) {
      extra<-matrix(0,length(forbidden),length(obj));extra[cbind(seq_along(forbidden),forbidden)]<-1
      mat<-rbind(mat,extra);dir<-c(dir,rep('==',length(forbidden)));rhs<-c(rhs,rep(0,length(forbidden)))
    }
    z<-lpSolve::lp(if(max)'max'else'min',obj,mat,ifelse(dir=='==','=',dir),rhs,
      all.int=TRUE,timeout=15L)
    solver_status<<-z$status
    if(z$status!=0L) stop(paste('matings exact allocation failed, lpSolve status',z$status))
    stopifnot(max(abs(z$solution-round(z$solution)))<1e-7)
    list(solution=round(z$solution),status=0L)
  }
  mating<-optisel_matings_dimension_safe(phen,cost,solver=solver_hook)
  result<-cbind(unname(index[mating$Dam]),unname(index[mating$Sire]))
  result<-result[rep(seq_len(nrow(result)),as.integer(mating$n)),,drop=FALSE]
  stopifnot(nrow(result)==crosses,all(allowed[result]),all(result[,1]!=result[,2]),
    identical(tabulate(plan[,1],n),tabulate(result[,1],n)),
    identical(tabulate(plan[,2],n),tabulate(result[,2],n)))
  list(plan=result,matings=mating,solver_status=solver_status,
    forbidden_edges=length(forbidden),mate_kinship=mean(K[result]))
}

# Greedy contribution repair; only starts from an optiSel-derived plan.
optisel_repair_gain <- function(plan,K,bv,allowed,B,max_parents=32L,max_steps=100L) {
  for(iter in seq_len(max_steps)) {
    met<-ocs_metrics(plan,K,bv)
    if(met$kinship<=B+1e-8) return(plan)
    best<-NULL;slots<-length(plan)
    for(slot in seq_len(slots)) {
      row<-(slot-1L)%%nrow(plan)+1L;col<-(slot-1L)%/%nrow(plan)+1L
      old<-plan[row,col];mate<-plan[row,3L-col];counts<-met$counts;counts[old]<-counts[old]-1L
      options<-which(if(col==1)allowed[,mate]else allowed[mate,])
      options<-options[counts[options]<ocs_cap_slots(nrow(plan))]
      if(sum(counts>0)>=max_parents)options<-options[counts[options]>0]
      c0<-counts/slots;Kc<-as.numeric(K%*%c0)
      kin<-sum(c0*Kc)+2*Kc[options]/slots+diag(K)[options]/slots^2
      reduction<-met$kinship-kin;ok<-reduction>1e-12
      options<-options[ok];kin<-kin[ok];reduction<-reduction[ok]
      if(!length(options))next
      loss<-(bv[old]-bv[options])/slots;score<-loss/pmin(reduction,met$kinship-B)
      j<-order(score,loss,kin,options)[1]
      if(is.null(best)||score[j]<best$score)best<-list(row=row,col=col,new=options[j],score=score[j])
    }
    if(is.null(best))break
    plan[best$row,best$col]<-best$new
  }
  # If greedy exchanges get stuck, decrease kinship locally from the same
  # optiSel-derived plan; never substitute a manual endpoint as the result.
  plan<-ocs_improve(plan,K,bv,allowed,max_parents,'diversity',sweeps=24L)
  if(ocs_metrics(plan,K,bv)$kinship>B+1e-8)stop('No optiSel-derived integer plan met the absolute kinship bound.')
  plan
}

optisel_contributions <- function(K,bv,bound=NULL,max_contribution=.10,N=1500) {
  phen<-data.frame(Indiv=rownames(K),Sex=NA_character_,BV=as.numeric(scale(bv)),isCandidate=TRUE)
  cand<-optiSel::candes(phen,N=N,Kin=K,quiet=TRUE)
  q<-cand$kinship$Kin
  stopifnot(max(abs(q$a))<1e-12)
  con<-list(ub=setNames(rep_len(max_contribution,nrow(K)),rownames(K)))
  if(!is.null(bound))con$ub.Kin<-bound+q$d
  fit<-optiSel::opticont(if(is.null(bound))'min.Kin'else'max.BV',cand,con,solver='slsqp',quiet=TRUE,maxeval=1500)
  if(!all(fit$info$valid))stop('optiSel continuous contribution constraint failed.')
  fit
}

optisel_matings_profile <- function(K,bv,allowed,crosses,max_parents,alpha,restarts=4L,sweeps=12L,seed=1L) {
  # Common feasible endpoint rule retained to compare profiles at identical
  # absolute bounds; manual gain-plan search is not used as an OCS result.
  had<-exists('.Random.seed',.GlobalEnv,inherits=FALSE)
  if(had)saved<-get('.Random.seed',.GlobalEnv)
  on.exit(if(had)assign('.Random.seed',saved,.GlobalEnv)else if(exists('.Random.seed',.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv))
  set.seed(seed)
  endpoints<-ocs_endpoints(K,bv,allowed,crosses,max_parents,restarts,sweeps)
  B<-endpoints$minK+alpha*(endpoints$maxK-endpoints$minK)
  maternal<-rowSums(allowed)>0;paternal<-colSums(allowed)>0
  role<-optisel_role_contributions(K,bv,maternal,paternal,B)
  fit<-role$fit;n<-length(bv);c<-numeric(n)
  sums<-rowsum(role$contributions,role$index,reorder=FALSE)
  c[as.integer(rownames(sums))]<-sums[,1]
  counts<-optisel_round_roles(role,crosses)
  rounded_plan<-NULL
  # Build an arbitrary plan with the requested role counts; matings() alone
  # decides its pairing, including hard zero bounds on forbidden edges.
  raw<-cbind(rep(counts$mothers,counts$mother_counts),rep(counts$fathers,counts$father_counts))
  allocated<-tryCatch(optisel_matings_exact(raw,K,allowed),error=function(e)e)
  recovery<-inherits(allocated,'error')
  if(!recovery) {
    plan<-allocated$plan;recovery<-ocs_metrics(plan,K,bv)$parents>max_parents || max(tabulate(plan,nbins=n))>ocs_cap_slots(crosses)
  }
  if(recovery) {
    z<-optisel_recover_transport(role,K,allowed,crosses,max_parents,bv,B,timeout=2L)
    plan<-z$plan
  }
  repaired<-ocs_metrics(plan,K,bv)$kinship>B+1e-8
  if(repaired)plan<-optisel_repair_gain(plan,K,bv,allowed,B,max_parents)
  plan<-ocs_improve(plan,K,bv,allowed,max_parents,'gain',B,sweeps=24L)
  final<-optisel_matings_exact(plan,K,allowed);plan<-final$plan
  met<-ocs_metrics(plan,K,bv)
  stopifnot(nrow(plan)==crosses,all(allowed[plan]),met$parents<=max_parents,met$kinship<=B+1e-8,max(met$contribution)<=.10+1e-10)
  list(plan=plan,metrics=met,bound=B,minK=endpoints$minK,maxK=endpoints$maxK,alpha=alpha,
    audit=list(continuous_gain=sum(c*bv),continuous_kinship=as.numeric(crossprod(c,K%*%c)),
      transport_recovery=recovery,kinship_repair=repaired,continuous_contributions=c,
      final_matings=final$matings,forbidden_edges=final$forbidden_edges,solver_status=final$solver_status,
      max_individual_contribution=.10,max_participations=ocs_cap_slots(crosses),
      implementation='optiSel::opticont + custom integer refinement + optiSel::matings (hard constraints)'))
}


# Different parental role eligibility (speed breeding) is explicit in the
# continuous fit: females and males each supply one half of contributions.
optisel_role_contributions <- function(K,bv,maternal,paternal,B) {
  m<-which(maternal);p<-which(paternal);n<-nrow(K)
  if(!length(m)||!length(p))stop('OCS has no eligible maternal or paternal role.')
  if(identical(m,p)) {
    fit<-optisel_contributions(K[m,m,drop=FALSE],bv[m],B)
    c<-fit$parent$oc[match(rownames(K)[m],fit$parent$Indiv)]
    return(list(fit=fit,contributions=c(c/2,c/2),index=c(m,p),mothers=seq_along(m),fathers=length(m)+seq_along(p)))
  }
  index<-c(m,p);ids<-c(paste0('M_',m),paste0('P_',p))
  Kr<-K[index,index,drop=FALSE];dimnames(Kr)<-list(ids,ids)
  scaled<-as.numeric(scale(bv))
  phen<-data.frame(Indiv=ids,Sex=c(rep('female',length(m)),rep('male',length(p))),BV=scaled[index],isCandidate=TRUE)
  cand<-optiSel::candes(phen,N=1500,Kin=Kr,quiet=TRUE);q<-cand$kinship$Kin
  stopifnot(max(abs(q$a))<1e-12)
  upper<-rep(.10,length(index));mi<-seq_along(m);pi<-length(m)+seq_along(p)
  for(attempt in 1:30) {
    fit<-optiSel::opticont('max.BV',cand,list(ub=setNames(upper,ids),ub.Kin=B+q$d),solver='slsqp',quiet=TRUE,maxeval=1500)
    if(!all(fit$info$valid))stop('Role-specific optiSel contribution constraints failed.')
    c<-pmax(0,fit$parent$oc[match(ids,fit$parent$Indiv)])
    combined<-numeric(n);sums<-rowsum(c,index,reorder=FALSE);combined[as.integer(rownames(sums))]<-sums[,1]
    stopifnot(all(is.finite(c)),abs(sum(c[mi])-.5)<1e-6,abs(sum(c[pi])-.5)<1e-6)
    dominant<-which(combined>.10+1e-7)
    if(!length(dominant))return(list(fit=fit,contributions=c,index=index,mothers=mi,fathers=pi,cap_attempts=attempt))
    for(i in dominant) {
      j<-which(index==i)
      upper[j]<-pmin(upper[j],.10*c[j]/combined[i])
    }
  }
  stop('Role-specific optiSel did not satisfy combined individual contribution caps.')
}

# Additional hard contribution cap and dimension-safe mate allocation.
source(file.path(getOption('almond.input_dir'),'ocs_cap10_support.R'),local=.GlobalEnv)
