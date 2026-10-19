# Focused allocation/integration tests only; no breeding experiment is launched.
# From the repository root: Rscript 02_ClonalBreeding/tests/check_ocs_native_offspring.R
root<-normalizePath('02_ClonalBreeding',winslash='/',mustWork=TRUE)
source(file.path(root,'ocs_discrete.R'));source(file.path(root,'ocs_optisel.R'))
# Small biology fixtures use a nonbinding bound; real captured inputs separately validate active bounds.
K<-diag(24);dimnames(K)<-list(as.character(1:24),as.character(1:24));bv<-seq_len(24)
rounding_warnings<-character()
for(kind in c('normal','speed'))for(Ne in c(100,50,25)) {
 allowed<-matrix(TRUE,24,24);diag(allowed)<-FALSE
 if(kind=='speed')allowed[9:24,]<-FALSE
 z<-withCallingHandlers(optisel_offspring_profile(K,bv,allowed,1500L,Ne,reference_kinship=.20,generation_interval=6),
   warning=function(w){rounding_warnings<<-c(rounding_warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 stopifnot(sum(z$family_sizes)==1500L,all(allowed[z$plan]),all(z$plan[,1]!=z$plan[,2]),
   max(z$metrics$contribution)<=.10+1e-10,z$audit$solver_status==0L,
   identical(z$audit$custom_contribution_refinement,FALSE),
   abs(sum(z$metrics$contribution)-1)<1e-10,
   abs(z$metrics$mate_kinship-weighted.mean(K[z$plan],z$family_sizes))<1e-10)
 n<-z$audit$offspring_counts
 stopifnot(sum(n$n)==3000L,z$audit$optimizer_variables==24L,
   abs(z$metrics$kinship-as.numeric(crossprod(z$metrics$contribution,K%*%z$metrics$contribution)))<1e-10)
 cat('PASS native profile:',kind,Ne,'parents',z$metrics$parents,'pairs',nrow(z$plan),'\n')
}
# Rounding departures are audited; no departure is required in a fixture.

# A deliberately cheapest forbidden edge must never be allocated.
phen<-data.frame(Indiv=c(paste0('M_',1:6),paste0('P_',7:12)),
 Sex=rep(c('female','male'),each=6),n=10L)
indices<-setNames(1:12,phen$Indiv)
K2<-diag(12);K2[1,7]<-K2[7,1]<- -100
A<-matrix(FALSE,12,12);A[1:6,7:12]<-TRUE;A[1,7]<-FALSE
z<-optisel_matings_offspring(phen,K2,A,indices)
stopifnot(sum(z$family_sizes)==60L,all(A[z$plan]),!any(z$plan[,1]==1 & z$plan[,2]==7))
# No private fallback when biological margins are impossible.
A[1,]<-FALSE
bad<-tryCatch(optisel_matings_offspring(phen,K2,A,indices),error=function(e)e)
stopifnot(inherits(bad,'error'),grepl('lpSolve status',conditionMessage(bad)))

# Verify variable family sizes at the crossing boundary using a mock cross engine.
# This checks IDs and offspring totals without simulating genetic generations.
if(!methods::isClass('AlmondCrossFixture'))
 methods::setClass('AlmondCrossFixture',slots=c(id='character',mother='character',father='character'))
methods::setMethod('[',signature(x='AlmondCrossFixture',i='ANY',j='missing',drop='missing'),
 function(x,i,j,...,drop)methods::new('AlmondCrossFixture',id=x@id[i],mother=x@mother[i],father=x@father[i]))
e<-new.env(parent=.GlobalEnv)
e$SP<-NULL;e$self_compatible_allele<-0;e$counter<-0L
e$nInd<-function(x)length(x@id)
e$ocs_compatibility<-function(pop,markers,self_compatible=NULL){a<-matrix(TRUE,length(pop@id),length(pop@id));diag(a)<-FALSE;a}
e$pullMarkerHaplo<-function(pop,...)matrix(0,length(pop@id),1)
e$makeCross<-function(pop,crossPlan,nProgeny,simParam=NULL){
 e$counter<-e$counter+nProgeny
 methods::new('AlmondCrossFixture',id=as.character(e$counter-seq_len(nProgeny)+1L),
 mother=rep(pop@id[crossPlan[1,1]],nProgeny),father=rep(pop@id[crossPlan[1,2]],nProgeny))}
e$mergePops<-function(xs)methods::new('AlmondCrossFixture',
 id=unlist(lapply(xs,function(x)x@id)),mother=unlist(lapply(xs,function(x)x@mother)),father=unlist(lapply(xs,function(x)x@father)))
sys.source(file.path(root,'compatible_crosses.R'),envir=e)
pop<-methods::new('AlmondCrossFixture',id=c('A','B','C'),mother=rep('',3),father=rep('',3))
state<-new.env();state$year<-10L;state$birth<-c(A=4L,B=6L,C=8L)
state$config<-list(minimum_mother_age=4L,minimum_pollen_age=2L)
previous<-getOption('almond.cross_state');options(almond.cross_state=state)
tryCatch({
 out<-e$randCrossGamSI(list(SIPos=1),pop,nCrosses=NA_integer_,nProgeny=75L,
   cross_plan=matrix(c('A','B','A','C'),ncol=2,byrow=TRUE),family_sizes=c(3L,7L))
 stopifnot(length(out@id)==10L,sum(out@father=='B')==3L,sum(out@father=='C')==7L,
   all(out@mother=='A'),identical(state$speed_cross_trace$family_sizes,c(3L,7L)),
   all(state$speed_cross_trace$mother_ages>=4),all(state$speed_cross_trace$father_ages>=2))
 bad<-tryCatch(e$randCrossGamSI(list(SIPos=1),pop,nCrosses=NA_integer_,nProgeny=75L,
   cross_plan=matrix(c('C','A'),ncol=2),family_sizes=10L),error=function(err)err)
 stopifnot(inherits(bad,'error'))
 # Non-OCS keeps its original scalar family size and crossing budget.
 out<-e$randCrossGamSI(list(SIPos=1),pop,nCrosses=2L,nProgeny=5L,
   cross_plan=matrix(c('A','B','A','C'),ncol=2,byrow=TRUE))
 stopifnot(length(out@id)==10L,identical(state$speed_cross_trace$family_sizes,c(5L,5L)))
},finally=options(almond.cross_state=previous))
cat('PASS: native profiles, 10% cap, forbidden edges, impossible margins, variable families and unchanged scalar families.\n')
