# Independent aggregate evaluations at the ends of PT, AT and elite trials.
# AlphaSimR generates each stage phenotype; individual stage records are cached
# so a genotype present in both a trial and the parental pool is evaluated once.

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

# Return the residual budget for a separate stage trial. The target total
# non-genetic variance includes one GxE environment, not an average of years.
stage_phenotype_budget <- function(state,stage) {
  endpoint <- if(stage=='Parents') 'ET3' else stage
  if(!endpoint %in% c('PT2','AT3','ET3')) stop('Not a completed phenotyping stage.')
  h2 <- unname(state$config$h2[[endpoint]])
  total <- state$reference_variance*(1/h2-1)
  residual <- total-state$reference_gxe_variance
  if(!is.finite(residual)||residual < -1e-8)
    stop('Stage target h2 is incompatible with the single-environment GxE variance.')
  list(endpoint=endpoint,total=total,gxe=state$reference_gxe_variance,
       residual=max(0,residual),h2=h2)
}

# Generate a completed-stage evaluation with native setPheno(). Intermediate
# years carry the last completed evaluation, or NA before completion of PT.
update_phenotypes <- function(pop,stage,state) {
  if(is.null(pop)) return(NULL)
  if(!identical(state$config$phenotyping_model,'independent_stage_setPheno_v1'))
    stop('Stage phenotyping requires a newly generated stage-model burn-in.')
  if(!is.environment(state$phenotype_records)) stop('Missing stage phenotype archive.')
  ids <- as.character(pop@id)
  endpoint <- if(stage=='Parents') 'ET3' else stage
  if(!endpoint %in% c('PT2','AT3','ET3')) {
    # Never retain a default phenotype generated when a new Pop was created.
    # Only completed trial records may be used by selection or prediction.
    for(i in seq_along(ids)) {
      records <- state$phenotype_records[[ids[i]]]
      available <- intersect(c('PT2','AT3','ET3'),names(records))
      if(!length(available)) pop@pheno[i,1] <- NA_real_ else {
        record <- records[[tail(available,1)]]
        pop@pheno[i,1] <- record$phenotype;pop@fixEff[i] <- record$year
      }
    }
    return(pop)
  }
  budget <- stage_phenotype_budget(state,endpoint)
  probability <- state$environment_p[[as.character(state$year)]]
  if(is.null(probability)||!is.finite(probability)) stop('Missing stage evaluation environment.')
  missing <- which(vapply(ids,function(id)is.null(state$phenotype_records[[id]][[endpoint]]),logical(1)))
  if(length(missing)) {
    evaluated <- AlphaSimR::setPheno(pop[missing],varE=budget$residual,reps=1,
      p=probability,fixEff=state$year,simParam=SP)
    for(j in seq_along(missing)) {
      i <- missing[j];id <- ids[i];records <- state$phenotype_records[[id]]
      if(is.null(records)) records <- list()
      records[[endpoint]] <- list(phenotype=as.numeric(evaluated@pheno[j,1]),
        year=state$year,stage=endpoint,p=probability,residual_variance=budget$residual,
        target_total_error=budget$total,target_h2=budget$h2)
      state$phenotype_records[[id]] <- records
    }
  }
  for(i in seq_along(ids)) {
    record <- state$phenotype_records[[ids[i]]][[endpoint]]
    pop@pheno[i,1] <- record$phenotype;pop@fixEff[i] <- record$year
  }
  pop
}

# Parents retain the last completed trial evaluation. At a trial endpoint they
# share the same cached record as their trial copy. Newly imported old founders
# receive an elite evaluation once; no independent annual re-phenotyping occurs.
update_parent_phenotypes <- function(state) {
  age <- state$year-state$birth[as.character(state$parents@id)]
  if(anyNA(age)) stop('Unknown parental birth year.')
  endpoint <- ifelse(age>=11L,'ET3',ifelse(age>=8L,'AT3',ifelse(age>=5L,'PT2','PT1')))
  parts <- lapply(unique(endpoint),function(stage) {
    update_phenotypes(state$parents[which(endpoint==stage)],stage,state)
  })
  state$parents <- mergePops(parts)
  invisible(state)
}
