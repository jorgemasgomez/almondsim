# Validate the active Ne-only configuration without simulating trajectories.
# Run from the repository root: Rscript 02_ClonalBreeding/tests/check_ne_configuration.R
root <- normalizePath('02_ClonalBreeding',winslash='/',mustWork=TRUE)
engine <- new.env(parent=.GlobalEnv)
sys.source(file.path(root,'redesigned_schemes.R'),envir=engine)
catalogue <- engine$redesign_all_scenarios(root)
stopifnot(nrow(catalogue)==18L,is.logical(catalogue$ocs_enabled),
          sum(catalogue$ocs_enabled)==6L)
targets <- c(Diversity=30,Balanced=10,Gain=5)
ocs <- catalogue[catalogue$ocs_enabled,,drop=FALSE]
stopifnot(identical(ocs$ocs_target_Ne,unname(targets[ocs$ocs_profile])),
          all(is.na(catalogue$ocs_target_Ne[!catalogue$ocs_enabled])))
for(i in seq_len(nrow(catalogue)))engine$validate_ocs_settings(as.list(catalogue[i,]))

default <- engine$redesign_config()
engine$validate_ocs_settings(default)
for(target in c(30,10,5))engine$validate_ocs_settings(list(ocs_enabled=TRUE,ocs_target_Ne=target))
invalid <- list(list(ocs_enabled=TRUE,ocs_target_Ne=NA_real_),
                list(ocs_enabled=TRUE,ocs_target_Ne=.5),
                list(ocs_enabled=TRUE,ocs_target_Ne=Inf),
                list(ocs_enabled=FALSE,ocs_target_Ne=30),
                list(ocs_target_Ne=30),list(ocs_enabled=1,ocs_target_Ne=30))
for(settings in invalid)stopifnot(inherits(try(engine$validate_ocs_settings(settings),silent=TRUE),'try-error'))
cat('PASS: the R catalogue and all 18 scheme directories agree; invalid OCS settings fail.\n')

# Replace prediction and optimization at the allocation boundary with spies.
# This tests dispatch and parental-use timing, not alternative genetic models.
if(!methods::isClass('NeDispatchFixture'))methods::setClass('NeDispatchFixture',
  slots=c(id='character',pheno='matrix',ebv='matrix'))
population <- function(ids)methods::new('NeDispatchFixture',id=ids,
  pheno=matrix(seq_along(ids)),ebv=matrix(seq_along(ids)))
spy <- new.env(parent=engine)
spy$redesign_predict <- function(pop,method,model)pop
spy$gv <- function(pop)seq_along(pop@id)
spy$called <- NULL
spy$redesign_ocs <- function(state,candidates,model) {
  spy$called <- 'OCS'
  list(parents=candidates,new=length(candidates@id),target_Ne=state$config$ocs_target_Ne)
}
spy$redesign_renew <- function(old,candidates,method,model,config) {
  spy$called <- 'replacement';list(parents=candidates,new=length(candidates@id))
}
allocate <- engine$redesign_allocate;environment(allocate)<-spy
for(i in seq_len(nrow(catalogue))) {
  settings <- as.list(catalogue[i,])
  state <- new.env();state$config <- settings;state$parents<-population(c('old1','old2','old3'))
  result <- allocate(state,population(c('new1','new2','new3')),settings$method,NULL,
                     settings$first_parental_use,41L)
  stopifnot(identical(spy$called,if(settings$ocs_enabled)'OCS' else 'replacement'),
            result$use_year==if(settings$first_parental_use>=6L)42L else 41L)
  if(settings$ocs_enabled)stopifnot(result$ocs$target_Ne==settings$ocs_target_Ne) else
    stopifnot(is.null(result$ocs))
}
cat('PASS: OCS versus replacement dispatch and current/next-year parent timing for all 18 schemes.\n')

# Reject incompatible source designs before modifying state or consuming RNG.
state <- new.env();state$config <- default;state$config$design_version <- 'unsupported_design'
set.seed(271L);seed <- .Random.seed
stopifnot(inherits(try(engine$redesign_year(state,'GS',6L,1L),silent=TRUE),'try-error'),
          identical(seed,.Random.seed),is.null(state$pop))
cat('PASS: unsupported source designs are rejected before annual state changes.\n')
