# Real AlphaSimR checks for the experimental, separate-stage phenotype model.
# No founder bank, long breeding trajectory or licensed solver is required.
library(AlphaSimR)
root <- normalizePath('02_ClonalBreeding',winslash='/',mustWork=TRUE)
source(file.path(root,'redesigned_schemes.R'))
set.seed(71021)
founders <- quickHaplo(nInd=120,nChr=1,segSites=200)
SP <- SimParam$new(founders)
SP$addTraitAG(nQtlPerChr=100,mean=100,var=150000,varGxE=150000)
SP$nThreads <- 1L
pop <- newPop(founders,simParam=SP)
state <- new.env(parent=emptyenv())
state$config <- breeding_config()
state$reference_variance <- 150000
state$reference_gxe_variance <- 150000
state$phenotype_records <- new.env(parent=emptyenv())
state$environment_p <- c(`1`=.3,`2`=.7,`5`=.6,`8`=.4)
state$year <- 1L
seed <- .Random.seed
pt1 <- update_phenotypes(pop,'PT1',state)
stopifnot(all(is.na(pt1@pheno)),identical(seed,.Random.seed),length(state$phenotype_records)==0)
state$year <- 2L
set.seed(456)
expected <- setPheno(pop,varE=450000,reps=1,p=.7,fixEff=2,simParam=SP)
set.seed(456)
pt <- update_phenotypes(pop,'PT2',state)
stopifnot(identical(pt@pheno,expected@pheno),identical(pt@fixEff,expected@fixEff))
seed <- .Random.seed
again <- update_phenotypes(pop,'PT2',state)
stopifnot(identical(pt@pheno,again@pheno),identical(seed,.Random.seed))
state$year <- 5L
carried <- update_phenotypes(pop,'AT2',state)
stopifnot(identical(pt@pheno,carried@pheno),all(carried@fixEff==2L),identical(seed,.Random.seed))
set.seed(789)
expected <- setPheno(pop,varE=75000,reps=1,p=.6,fixEff=5,simParam=SP)
set.seed(789)
at <- update_phenotypes(pop,'AT3',state)
stopifnot(identical(at@pheno,expected@pheno),!identical(pt@pheno,at@pheno))
stopifnot(all(vapply(as.list(state$phenotype_records),function(x)
  identical(names(x),c('PT2','AT3'))&&x$PT2$year==2L&&x$AT3$year==5L,logical(1))))
state$parents <- pop
state$birth <- setNames(rep(-3L,nInd(pop)),as.character(pop@id))
seed <- .Random.seed
update_parent_phenotypes(state)
stopifnot(identical(state$parents@pheno,at@pheno),identical(seed,.Random.seed))
state$year <- 8L
elite <- update_phenotypes(pop,'ET3',state)
stopifnot(all(is.finite(elite@pheno)),stage_phenotype_budget(state,'ET3')$residual==0)
calibration <- calculate_pt_at_residual_covariance(state)
stopifnot(identical(unname(calibration$absolute_working_matrix),diag(c(600000,225000))),
          calibration$residual_correlation==0,
          identical(unname(calibration$relative_matrix),diag(c(1,.375))))
# Completed records are archived before selection; intermediate years create
# no apparent extra PT/AT observations and no look-ahead data.
state$stage_history <- list();state$pop <- list(PT1=pt1,PT2=pt,AT1=carried,AT3=at,ET3=elite)
archive_evaluated_cohorts(state)
stopifnot(setequal(names(state$stage_history),c('PT2','AT3','ET3')))
cat('PASS: native setPheno equivalence, absent pre-PT records, endpoint caching, parent/trial identity, separate PT/AT records, zero cross-stage error covariance, and completed-cohort archives.\n')
