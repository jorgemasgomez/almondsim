# Isolated prototype: direct bound, single opticont call, default solver.
# No simulation/configuration changes. Ne conversion requires calibrated reference and generation interval.
# Native optiSel contribution -> offspring count -> mating workflow.
# This adapter defines almond biology and validates results; it does not repair,
# locally refine or substitute contribution vectors returned by optiSel.
if(!requireNamespace('optiSel',quietly=TRUE)) {
 lib<-Sys.getenv('ALMONDSIM_OCS_R_LIBRARY','C:/Users/franc/Documents/Codex/2026-10-04/es/work/r-library')
 if(dir.exists(lib)).libPaths(c(lib,.libPaths()))
}
if(!requireNamespace('optiSel',quietly=TRUE)||!requireNamespace('lpSolve',quietly=TRUE))
 stop('OCS requires optiSel and lpSolve.')

optisel_offspring_context <- function(K,bv,allowed,N,max_contribution=.10,si_groups=NULL,individual_caps=NULL) {
 stopifnot(is.matrix(K),nrow(K)==ncol(K),length(bv)==nrow(K),identical(dim(allowed),dim(K)),
 !anyNA(allowed),all(is.finite(K)),all(is.finite(bv)),N>=1,N==floor(N),max_contribution>0,max_contribution<=1,
 max(abs(K-t(K)))<1e-8)
 diag(allowed)<-FALSE
 eligible<-rowSums(allowed)>0|colSums(allowed)>0
 if(is.null(individual_caps))individual_caps<-rep(max_contribution,length(bv))
 stopifnot(length(individual_caps)==length(bv),all(is.finite(individual_caps)),all(individual_caps>0),all(individual_caps<=max_contribution))
 if(!any(eligible))stop('No eligible OCS candidates.')
 ids<-paste0('G_',seq_along(bv));dimnames(K)<-list(ids,ids)
 phen<-data.frame(Indiv=ids,Sex=NA_character_,BV=if(sd(bv)>0)as.numeric(scale(bv))else 0,
 isCandidate=eligible,OnlyPollen=as.numeric(rowSums(allowed)==0),OnlyMother=as.numeric(colSums(allowed)==0))

 group_bounds<-list()
 if(!is.null(si_groups)) {
  stopifnot(length(si_groups)==length(bv),!anyNA(si_groups))
  gs<-split(seq_along(bv),si_groups)
  for(j in seq_along(gs)) if(length(gs[[j]])*max_contribution>.5+1e-10) {
   ii<-gs[[j]];stopifnot(!any(allowed[ii,ii,drop=FALSE]))
   nm<-paste0('SIgroup',j);phen[[nm]]<-as.numeric(seq_along(bv)%in%ii);group_bounds[[paste0('ub.',nm)]]<-.5
  }
 }
 # One optimizer variable per genotype. Native trait bounds restrict role-only groups.
 cand<-optiSel::candes(phen,N=N,Kin=K,quiet=TRUE)
 q<-cand$kinship$Kin
 stopifnot(max(abs(q$a))<1e-10)
 list(K=K,bv=bv,allowed=allowed,N=as.integer(N),cap=individual_caps,
 si_groups=si_groups,cand=cand,con=c(list(ub=setNames(ifelse(eligible,individual_caps,0),ids),ub.OnlyPollen=.5,ub.OnlyMother=.5),group_bounds),
 drift=q$d,index=seq_along(bv),ids=ids)
}

optisel_native_fit <- function(context,method,bound=NULL) {
 con<-context$con
 if(!is.null(bound))con$ub.Kin<-bound+context$drift
 fit<-optiSel::opticont(method,context$cand,con,solver='default',quiet=TRUE)
 if(!is.null(getOption('almond.decision_capture_dir')))saveRDS(list(info=fit$info,fit=fit,bound=bound,caps=context$cap),file.path(getOption('almond.decision_capture_dir'),'current_native_fit_before_validation.rds'))
 if(!all(fit$info$valid))stop('Native optiSel continuous constraints failed.')
 oc<-fit$parent$oc[match(context$ids,fit$parent$Indiv)]
 if(anyNA(oc)||any(!is.finite(oc))||any(oc< -1e-7))stop('Invalid optiSel contribution vector.')
 combined<-numeric(nrow(context$K))
 for(j in seq_along(oc))combined[context$index[j]]<-combined[context$index[j]]+oc[j]
 if(!is.null(context$si_groups)&&max(tapply(combined,context$si_groups,sum))>.5+1e-7)stop('Continuous SI group contribution exceeds50%.')
 if(abs(sum(combined)-1)>1e-6||any(combined>context$cap+1e-7))
   stop('Native optiSel combined individual contribution cap failed.')
 kin<-as.numeric(crossprod(combined,context$K%*%combined))
 continuous_excess<-if(is.null(bound))0 else max(0,kin-bound)
 continuous_relative<-if(is.null(bound)||continuous_excess==0)0 else if(abs(bound)>1e-12)continuous_excess/abs(bound)else Inf
 if(!is.null(bound)&&continuous_excess>.10*abs(bound)+1e-7)
  stop(sprintf('Continuous coancestry exceeds authorized 10%% relative allowance: %.6g%%.',100*continuous_relative))
 if(continuous_excess>1e-7)warning(sprintf('Continuous coancestry excess %.6g%% accepted under user-authorized 10%% policy.',100*continuous_relative),call.=FALSE)
 list(fit=fit,contribution=combined,kinship=kin,gain=sum(combined*context$bv),continuous_excess=continuous_excess,continuous_relative=continuous_relative)
}

optisel_matings_offspring <- function(phen,K,allowed,index_by_id) {
 if(!is.null(getOption('almond.decision_capture_dir')))saveRDS(list(year=getOption('almond.timing_year'),phen=phen,K=K,allowed=allowed,index_by_id=index_by_id),file.path(getOption('almond.decision_capture_dir'),'current_mating_input.rds'))
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

# Allocate integer TOTAL participations jointly to directed mating edges.
# This is an explicit lpSolve adaptation, not optiSel::matings(). It preserves
# native noffspring counts exactly and does not repair/refine contributions.
optisel_joint_roles <- function(counts,K,allowed,N) {
 if(!is.null(getOption('almond.decision_capture_dir')))saveRDS(list(year=getOption('almond.timing_year'),counts=counts,K=K,allowed=allowed,N=N),file.path(getOption('almond.decision_capture_dir'),'current_joint_mating_input.rds'))
 active<-which(counts>0);edges<-which(allowed & outer(counts>0,counts>0,'&'),arr.ind=TRUE)
 if(!nrow(edges))stop('No compatible edges among contributing genotypes.')
 mat<-matrix(0,length(active)+1L,nrow(edges))
 for(j in seq_along(active))mat[j,]<-as.integer(edges[,1]==active[j])+as.integer(edges[,2]==active[j])
 mat[nrow(mat),]<-1
 z<-lpSolve::lp('min',K[edges],mat,rep('=',nrow(mat)),c(counts[active],N),all.int=TRUE,timeout=30L)
 if(z$status!=0L)stop('Joint role/mating allocation infeasible or unfinished: lpSolve status ',z$status)
 stopifnot(max(abs(z$solution-round(z$solution)))<1e-7)
 used<-which(z$solution>.5)
 list(plan=edges[used,,drop=FALSE],family_sizes=as.integer(round(z$solution[used])),solver_status=z$status,
 matings=NULL,forbidden_edges=sum(!allowed),backend='lpSolve joint roles and matings')
}

optisel_offspring_profile <- function(K,bv,allowed,N,target_Ne,max_contribution=.10,reference_kinship,generation_interval,si_groups=NULL,age=NULL,adult_cap=max_contribution) {
 stopifnot(length(target_Ne)==1L,is.finite(target_Ne),target_Ne>.5,
 is.finite(reference_kinship),reference_kinship<1,is.finite(generation_interval),generation_interval>0)
 stopifnot(length(age)==length(bv),all(is.finite(age)),adult_cap>0,adult_cap<=max_contribution)
 adult<-age>=4;individual_caps<-ifelse(adult,adult_cap,max_contribution)
 minimum_adults<-if(adult_cap<max_contribution)ceiling((if(all(adult))1 else .5)/adult_cap)else 0L
 if(sum(adult & rowSums(allowed)>0)<minimum_adults)stop(sprintf('Insufficient eligible adults for3.2%% cap: available%d, required%d.',sum(adult & rowSums(allowed)>0),minimum_adults))
 # Tutorial-style annual bound. Ne is a control parameter; founder-centred K
 # and standing-parent age proxy do not imply a calibrated realized effective size.
 B<-1-(1-reference_kinship)*(1-1/(2*target_Ne))^(1/generation_interval)
 started<-proc.time()
 context<-optisel_offspring_context(K,bv,allowed,N,max_contribution,si_groups,individual_caps)
 context_seconds<-unname((proc.time()-started)['elapsed']);fit_started<-proc.time()
 optimum<-optisel_native_fit(context,'max.BV',B)
 if(!is.null(getOption('almond.decision_capture_dir')))saveRDS(list(year=getOption('almond.timing_year'),optimum=optimum,si_groups=si_groups),file.path(getOption('almond.decision_capture_dir'),'current_continuous_fit.rds'))
 optimization_seconds<-unname((proc.time()-fit_started)['elapsed']);allocation_started<-proc.time()
 c<-optimum$contribution;ids<-context$ids
 only_pollen<-rowSums(context$allowed)==0;only_mother<-colSums(context$allowed)==0
 if(!any(only_pollen|only_mother)) {
  # Equal roles are feasible for the continuous age constraints in GS6.
  phen<-data.frame(Indiv=c(paste0('M_',seq_along(c)),paste0('P_',seq_along(c))),
   Sex=rep(c('female','male'),each=length(c)),oc=rep(c/2,2))
  count<-optiSel::noffspring(phen,N=N,random=FALSE);phen$n<-count$nOff[match(phen$Indiv,count$Indiv)]
  index<-setNames(rep(seq_along(c),2),phen$Indiv)
  rounded<-as.numeric(rowsum(phen$n,rep(seq_along(c),2),reorder=FALSE))/(2*N)
  if(any(rounded>individual_caps+1e-10))stop('Native rounding exceeded individual cap; no repair applied.')
  final<-optisel_matings_offspring(phen,K,context$allowed,index)
  mode<-'native equal-role noffspring and matings'
 } else if(any(only_pollen)&&abs(sum(c[only_pollen])-.5)<1e-7) {
  # Saturated pollen-only group: all other contributions must be maternal.
  phen<-data.frame(Indiv=ids,Sex=ifelse(only_pollen,'male','female'),oc=c)
  count<-optiSel::noffspring(phen,N=N,random=FALSE);phen$n<-count$nOff[match(ids,count$Indiv)]
  rounded<-phen$n/(2*N)
  if(any(rounded>individual_caps+1e-10))stop('Native rounding exceeded individual cap; no repair applied.')
  final<-optisel_matings_offspring(phen,K,context$allowed,setNames(seq_along(c),ids))
  mode<-'native saturated pollen-only group noffspring and matings'
 } else {
  # A bookkeeping group lets native noffspring allocate 2N total participations.
  # Biological roles remain free until the directed integer allocation.
  phen<-data.frame(Indiv=ids,Sex='male',oc=c/2)
  count<-optiSel::noffspring(phen,N=2L*N,random=FALSE);phen$n<-count$nOff[match(ids,count$Indiv)]
  rounded<-phen$n/(2*N)
  if(any(rounded>individual_caps+1e-10))stop('Native rounding exceeded individual cap; no repair applied.')
  final<-optisel_joint_roles(phen$n,K,context$allowed,N)
  mode<-'native total noffspring then joint lpSolve roles/matings'
 }
 if(anyNA(phen$n)||any(phen$n<0)||any(phen$n!=floor(phen$n))||sum(phen$n)!=2*N)
  stop('Native offspring counts did not preserve 2N integer participations.')
 met<-ocs_metrics(final$plan,K,bv,weights=final$family_sizes)
 stopifnot(final$solver_status==0L,sum(final$family_sizes)==N,all(context$allowed[final$plan]),
 all(final$plan[,1]!=final$plan[,2]),all(met$contribution<=individual_caps+1e-10),sum(adult & met$counts>0)>=minimum_adults,
 max(abs(met$contribution-rounded))<1e-10)
 excess<-max(0,met$kinship-B);relative<-if(B>0)excess/B else if(excess==0)0 else Inf
 if(excess>.10*abs(B)+1e-7)stop('Rounded coancestry exceeds authorized 10% relative allowance.')
 if(excess>1e-8)warning(sprintf('Native offspring rounding coancestry exceeds target by %.12g (%.6g%%); counts preserved.',excess,100*relative),call.=FALSE)
 list(plan=final$plan,family_sizes=final$family_sizes,metrics=met,bound=B,minK=NA_real_,maxK=NA_real_,alpha=NA_real_,
 audit=list(implementation='optiSel SI groups adultcap3.2percent v12',individual_caps=individual_caps,adult_age_threshold=4,adult_cap=adult_cap,minimum_adults=minimum_adults,realized_adult_contributors=sum(adult & met$counts>0),si_groups=si_groups,SI_group_bound=.5,max_continuous_SI_group=if(is.null(si_groups))NA_real_ else max(tapply(c,si_groups,sum)),allocation_mode=mode,
 coancestry_policy='relative_excess_max_10pct_20261018',max_relative_coancestry_excess=.10,continuous_kinship_excess=optimum$continuous_excess,continuous_kinship_excess_relative=optimum$continuous_relative,optimizer_variables=length(c),context_seconds=context_seconds,optimization_seconds=optimization_seconds,allocation_seconds=unname((proc.time()-allocation_started)['elapsed']),target_Ne=target_Ne,reference_kinship=reference_kinship,generation_interval_proxy=generation_interval,
 continuous_gain=optimum$gain,continuous_kinship=optimum$kinship,continuous_contributions=c,
 continuous_fit=optimum$fit,rounded_contributions=rounded,offspring_counts=phen,
 final_matings=final$matings,solver_status=final$solver_status,continuous_bound=B,
 max_individual_contribution=max_contribution,offspring=N,rounded_bound_satisfied=excess<=1e-8,
 rounded_kinship_excess=excess,rounded_kinship_excess_relative=relative,
 max_parents=NA_integer_,custom_contribution_refinement=FALSE,custom_rounding_allowance=FALSE))
}
