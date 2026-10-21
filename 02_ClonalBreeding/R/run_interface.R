# Human-facing names map to the existing scenario IDs; scientific identifiers
# and result columns are unchanged. E denotes early selection, not a new method.
scheme_aliases <- function() {
  c(PH9='Pheno_Y9', PED9='Pedigree_Y9', GS9='GS_Y9',
    PH6='Pheno_Y6', PED6='Pedigree_Y6', GS6='GS_Y6',
    PED4='Pedigree_Y4', GS4='GS_Y4', SB2='GS_SpeedBreeding_Y2',
    'PED6-E'='Pedigree_Y6_Early1500', 'GS6-E'='GS_Y6_Early1500',
    'SB2-E'='GS_SpeedBreeding_Y2_Early1500',
    'GS6-OD'='GS_Y6_Early1500_OCS_Diversity',
    'GS6-OB'='GS_Y6_Early1500_OCS_Balanced',
    'GS6-OG'='GS_Y6_Early1500_OCS_Gain',
    'SB2-OD'='GS_SpeedBreeding_Y2_Early1500_OCS_Diversity',
    'SB2-OB'='GS_SpeedBreeding_Y2_Early1500_OCS_Balanced',
    'SB2-OG'='GS_SpeedBreeding_Y2_Early1500_OCS_Gain')
}

# Select schemes without changing catalogue order (and therefore worker order).
select_breeding_schemes <- function(root, schemes) {
  catalogue <- all_breeding_scenarios(root)
  aliases <- scheme_aliases()
  ids <- vapply(schemes, function(id) {
    if (id %in% names(aliases)) unname(aliases[[id]]) else id
  }, character(1))
  unknown <- setdiff(ids, catalogue$scenario_id)
  if (length(unknown)) stop('Unknown breeding scheme(s): ', paste(unknown, collapse=', '))
  catalogue[catalogue$scenario_id %in% ids, , drop=FALSE]
}

# The same engine serves one scheme or a selected experiment. This interface
# adds validation only; it does not implement a second simulation algorithm.
run_breeding_experiment <- function(root, schemes='all', config=breeding_config()) {
  scenarios <- if (identical(schemes,'all')) all_breeding_scenarios(root) else
    select_breeding_schemes(root, schemes)
  if (!nrow(scenarios)) stop('Select at least one breeding scheme.')
  coordinate_breeding_experiment(root, config, scenarios)
}

run_breeding_scheme <- function(root, scheme, config=breeding_config()) {
  if (length(scheme)!=1L) stop('run_breeding_scheme requires exactly one scheme.')
  run_breeding_experiment(root, scheme, config)
}

# CLI configuration is explicit. Help and validation do not create a run.
pipeline_cli <- function(root, single=FALSE, args=commandArgs(TRUE)) {
  usage <- c(
    'Run from the repository root after setting ALMONDSIM_FOUNDER_BANK and ALMONDSIM_DATA_DIR.',
    if(single) 'Rscript 02_ClonalBreeding/run_scheme.R --scheme=GS6-E --replicates=1 --workers=1 --threads=6' else
      'Rscript 02_ClonalBreeding/run_experiment.R --schemes=GS6-E,SB2-E --replicates=10 --workers=4 --threads=6',
    if(single) 'Use --list to print aliases without starting simulations.' else
      'Use --schemes=all for all 18 schemes. Use --list to print aliases without starting simulations.',
    'Each execution starts a NEW experiment. These launchers do not resume or retry old runs.'
  )
  if (!length(args) || '--help' %in% args) {cat(paste(usage,collapse='\n'),'\n');return(invisible(NULL))}
  if (identical(args,'--list')) {print(scheme_aliases());return(invisible(NULL))}
  if (any(!grepl('^--[a-z]+=.+$',args))) stop('Use --option=value arguments; see --help.')
  keys <- sub('^--([^=]+)=.*$','\\1',args)
  values <- sub('^--[^=]+=','',args)
  allowed <- c(if(single)'scheme' else 'schemes','replicates','workers','threads')
  if (anyDuplicated(keys) || any(!keys %in% allowed)) stop('Unknown or duplicated CLI option; see --help.')
  options <- setNames(as.list(values),keys)
  selection <- options[[if(single)'scheme' else 'schemes']]
  if (is.null(selection)) stop('Specify the scheme(s) explicitly; see --help.')
  if(single && (selection=='all'||grepl(',',selection))) stop('Select exactly one scheme.')
  config <- breeding_config()
  for (key in intersect(keys,c('replicates','workers','threads'))) {
    value <- suppressWarnings(as.numeric(options[[key]]))
    if (!is.finite(value)||value<1||value!=floor(value)||value>.Machine$integer.max)
      stop(key,' must be a positive integer.')
    config[[if(key=='replicates')'reps' else key]] <- as.integer(value)
  }
  run_breeding_experiment(root, strsplit(selection,',',fixed=TRUE)[[1]], config)
}
