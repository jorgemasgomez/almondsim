# Validate and import the founder bank, initialize shared state, and fill the pipeline before burn-in.
# Function names describe each operation; all launchers use this shared implementation.

# Locate a completed 150-diploid, eight-chromosome bank and reject incomplete or incompatible inputs.
find_founder_bank <- function(bank_root) {
  paths <- list.files(bank_root,pattern='^run_',full.names=TRUE)
  paths <- paths[dir.exists(paths)]
  if(file.exists(file.path(bank_root,'manifest.json'))) paths <- c(bank_root,paths)
  if(length(paths)) paths <- paths[order(file.info(paths)$mtime,decreasing=TRUE)]
  for(path in paths) {
    manifest_file <- file.path(path,'manifest.json')
    if(!file.exists(manifest_file)) next
    manifest <- jsonlite::fromJSON(manifest_file)
    if(!identical(manifest$status,'complete')) next
    if(manifest$configuration$n_diploid_individuals!=150L ||
       manifest$configuration$mutation_rate!=1e-8 || nrow(manifest$chromosomes)!=8L ||
       !setequal(as.character(manifest$chromosomes$chromosome),as.character(1:8)) ||
       anyNA(manifest$chromosomes$exported_sites) ||
       any(manifest$chromosomes$exported_sites<=0)) next
    files <- c(paste0('Chr_',1:8,'_position.txt'),paste0('Chr_',1:8,'_gmatrix.txt'))
    if(all(file.exists(file.path(path,files)))) return(path)
  }
  stop('No complete bank of 150 diploid founders with eight chromosomes and mutation_rate=1e-8 is available yet. Wait for bank150 generation to finish.')
}

# Validate the required founder parent count after pipeline initialization; do not invent additional parents.
validate_founder_parent_pool <- function(state) {
 if(nInd(state$parents)!=state$config$parents) stop('The revised design requires exactly 32 sampled founders, not a legacy pool.')
 invisible(state)
}

# Initialize AlphaSimR globals and state, then perform fill cycles without parental renewal.
initialize_breeding_pipeline <- function(root,config) {
  # Reuse validated founder import and S-locus assignment from the original program.
  source(file.path(root,'00_Burn_in','GlobalParameters.R'),local=.GlobalEnv)
  source(file.path(root,'compatible_crosses.R'),local=.GlobalEnv)
  source(file.path(root,'ocs_discrete.R'),local=.GlobalEnv)
  source(file.path(root,'optimize_contributions.R'),local=.GlobalEnv)
  if(!identical(config$founder_source,'bank150')) stop('The revised design requires the bank150/60000 founder source.')
  bank <- find_founder_bank(config$founder_bank)
  founder_count <- config$parents
  assign('nParents',founder_count,envir=.GlobalEnv)
  options(almond.haplotype_dir=bank,almond.founder_count=founder_count)
  source(file.path(root,'CreateFounders_redesigned.R'),local=.GlobalEnv)
  stopifnot(nInd(Parents)==founder_count)
  state <- new.env(parent=emptyenv())
  state$config <- config
  state$founder_count <- founder_count
  state$founder_indices <- founder_indices
  state$ocs_frequency <- ocs_fixed_frequency(Parents)
  state$ocs_trace <- list()
  state$reference_variance <- initVarG
  state$reference_gxe_variance <- initVarGE
  state$annual_evaluations <- new.env(parent=emptyenv())
  environment_seed <- if(is.null(config$environment_seed)) config$seed+900000L else config$environment_seed
  state$environment_p <- create_environment_schedule(config,environment_seed)
  state$errors <- new.env(parent=emptyenv())
  state$parents <- Parents
  state$birth <- setNames(rep(-30L,nInd(Parents)),as.character(Parents@id))
  state$pop <- setNames(vector('list',length(config$stages)),config$stages)
  state$stage_history <- list()
  state$training_stage <- NULL
  state$pending <- NULL
  state$prediction_model <- NULL
  state$model_method <- NULL
  state$model_training_stage <- NULL
  state$year <- -config$fill_cycles
  state$parents <- update_phenotypes(state$parents,'Parents',state)
  for(i in seq_len(config$fill_cycles)) {
    state$year <- i-config$fill_cycles
    advance_breeding_year(state,'Pheno',9L,0L,record=FALSE,renew=FALSE)
  }
  validate_founder_parent_pool(state)
  state
}

