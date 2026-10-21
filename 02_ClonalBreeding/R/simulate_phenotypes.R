# Generate aggregate yield evaluations and their correlated errors. Reference variances remain fixed after selection.
# Function names describe each operation; all launchers use this shared implementation.

# Update a correlated aggregate residual while preserving the requested marginal variance.
update_correlated_error <- function(previous,old_variance,new_variance) {
  if(old_variance<=1e-12) return(rnorm(length(previous),sd=sqrt(new_variance)))
  covariance <- min(old_variance,new_variance)
  ratio <- covariance/old_variance
  ratio*previous+rnorm(length(previous),sd=sqrt(max(0,new_variance-covariance^2/old_variance)))
}

# Create the shared annual environment calendar and restore the caller's RNG state.
create_environment_schedule <- function(config,seed) {
  had_seed <- exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) old_seed <- get('.Random.seed',envir=.GlobalEnv)
  on.exit(if(had_seed) assign('.Random.seed',old_seed,envir=.GlobalEnv) else
    if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv))
  set.seed(seed)
  years <- (-config$fill_cycles):(config$burnin+config$future)
  setNames(runif(length(years)),as.character(years))
}

# Update one population's aggregate phenotype using its evaluation history and stage precision.
update_phenotypes <- function(pop,stage,state) {
  if(is.null(pop)) return(NULL)
  if(is.null(state$environment_p) || is.null(state$annual_evaluations))
    stop('Missing annual environment state. Generate a new burn-in.')
  probability <- state$environment_p[[as.character(state$year)]]
  if(is.null(probability) || !is.finite(probability)) stop('Annual environment is missing for this year.')
  trait <- SP$traits[[1]]
  if(!.hasSlot(trait,'envVar') || !length(pop@gxe)) stop('Annual evaluation requires the original GxE trait.')
  annual_z <- qnorm(probability,sd=sqrt(trait@envVar))
  residuals <- environment <- numeric(nInd(pop))
  ids <- as.character(pop@id)
  for(i in seq_along(ids)) {
    history <- state$annual_evaluations[[ids[i]]]
    if(is.null(history)) history <- list(years=integer(),sum_z=0,n=0L)
    if(!state$year %in% history$years) {
      history$years <- c(history$years,state$year)
      history$sum_z <- history$sum_z+annual_z;history$n <- history$n+1L
      state$annual_evaluations[[ids[i]]] <- history
    }
    environment[i] <- pop@gxe[[1]][i]*history$sum_z/history$n
    # Budget includes averaged GxE: do not add it on top of the full target error.
    total_error <- state$reference_variance*(1/state$config$h2[[stage]]-1)
    variance <- total_error-state$reference_gxe_variance/history$n
    if(variance < -1e-8) stop('Target aggregate h2 is incompatible with annual GxE variance.')
    variance <- max(0,variance)
    prior <- state$errors[[ids[i]]]
    if(is.null(prior)) error <- rnorm(1,sd=sqrt(variance)) else
      error <- update_correlated_error(prior$error,prior$variance,variance)
    state$errors[[ids[i]]] <- list(error=error,variance=variance)
    residuals[i] <- error
  }
  pop@pheno <- matrix(as.numeric(gv(pop))+environment+residuals,ncol=1)
  pop@fixEff <- rep(as.integer(state$year),nInd(pop))
  pop
}

# Evaluate retained parents according to their ages; juvenile parents have no yield phenotype.
update_parent_phenotypes <- function(state) {
 age <- state$year-state$birth[as.character(state$parents@id)]
 if(anyNA(age)) stop('Unknown parental birth year.')
 parts <- list();juvenile <- which(age<4L)
 if(length(juvenile)) {
  pop <- state$parents[juvenile];pop@pheno[,] <- NA_real_;parts[[1L]] <- pop
 }
 stages <- setdiff(names(state$config$h2),'Parents')
 counts <- pmin(length(stages),pmax(1L,age-3L))
 for(count in sort(unique(counts[age>=4L]))) {
  pop <- state$parents[which(age>=4L & counts==count)]
  parts[[length(parts)+1L]] <- update_phenotypes(pop,stages[count],state)
 }
 state$parents <- mergePops(parts)
 invisible(state)
}

