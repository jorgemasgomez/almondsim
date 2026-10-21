# Define shared parameters and the 18 supported schemes. Loading this file does not run simulations.
# Function names describe each operation; all launchers use this shared implementation.

# Return defaults; environment variables supply external data and founder locations.
breeding_config <- function() {
 list(design_version='PT_AT_ET_optiSel_adultcap32_v12',reps=10L,workers=4L,threads=6L,
  seed=20261005L,asreml_workspace='16gb',burnin=40L,future=40L,fill_cycles=14L,
  parents=32L,crosses=20L,progeny=25L,act=25L,ect=8L,entrants=8L,
  training_cohorts=3L,gs_solver='AlphaSimR::RRBLUP',early_selection=FALSE,seedlings_keep=500L,
  ocs_enabled=FALSE,ocs_target_Ne=NA_real_,ocs_max_contribution=.10,ocs_backend='optiSel_adultcap32_v12',
  speed_breeding=FALSE,minimum_mother_age=4L,minimum_pollen_age=4L,
  data_dir=Sys.getenv('ALMONDSIM_DATA_DIR',unset=''),founder_source='bank150',
  founder_bank=Sys.getenv('ALMONDSIM_FOUNDER_BANK',unset=
    'C:/Users/franc/Desktop/Proyectos/Almond_sim/Haplotipos_2026-10-05_60000'),
  stages=c('F1','Seedlings','Juvenile2','Juvenile3','PT1','PT2',
           'AT1','AT2','AT3','ET1','ET2','ET3'),
  h2=c(PT1=.2,PT2=.2,AT1=.4,AT2=.4,AT3=.4,ET1=.5,ET2=.5,ET3=.5,Parents=.5),
  parent_stage=c(`9`='AT3',`6`='PT2',`4`='PT1',`2`='Juvenile2'),
  training_stage=c(`9`='PT2',`6`='PT2',`4`='PT2',`2`='PT2'),
  training_cohort_stage=c(`9`='AT3',`6`='PT2',`4`='PT2',`2`='PT2'))
}

# An explicit switch activates OCS; the target Ne alone controls its bound.
# Reject incomplete settings instead of silently treating them as non-OCS.
validate_ocs_settings <- function(config) {
 if(!is.logical(config$ocs_enabled)||length(config$ocs_enabled)!=1L||is.na(config$ocs_enabled))
   stop('ocs_enabled must be TRUE or FALSE.')
 target <- config$ocs_target_Ne
 if(!is.numeric(target)||length(target)!=1L)
   stop('ocs_target_Ne must be one numeric value (NA when OCS is disabled).')
 if(config$ocs_enabled && (!is.finite(target)||target<=.5))
   stop('Enabled OCS requires a finite ocs_target_Ne greater than 0.5.')
 if(!config$ocs_enabled && !is.na(target))
   stop('Disabled OCS requires ocs_target_Ne = NA_real_.')
 invisible(TRUE)
}

# Build the scientific catalogue and optionally validate per-scheme GlobalParameters.R files.
breeding_scenarios <- function(root=NULL) {
 make <- function(method,first,early=FALSE,renewal=25,ocs_enabled=FALSE,profile=NA_character_) {
  speed <- first==2L
  suffix <- if(ocs_enabled) paste0('_OCS_',profile) else if(renewal!=25) paste0('_R',renewal) else ''
  stem <- paste0(method,if(speed) '_SpeedBreeding' else '', '_Y',first,
                 if(early) '_Early1500' else '',suffix)
  data.frame(method=method,first_parental_use=as.integer(first),
   phenotyping_years=if(first==9L) 5L else if(first==6L) 2L else 0L,
   parent_selection_stage=breeding_config()$parent_stage[[as.character(first)]],
   training_stage=breeding_config()$training_stage[[as.character(first)]],
   training_cohort_stage=breeding_config()$training_cohort_stage[[as.character(first)]],
   selection_timing=if(first>=6L) 'end_of_previous_year' else 'before_own_phenotyping',
   early_selection=early,progeny_per_cross=if(early) 75L else 25L,
   seedlings_keep=500L,new_parents_per_year=if(!ocs_enabled) as.integer(32*renewal/100) else NA_integer_,
   renewal_percent=if(!ocs_enabled) renewal else NA_real_,ocs_enabled=ocs_enabled,
   ocs_profile=profile,ocs_target_Ne=if(!ocs_enabled) NA_real_ else unname(c(Diversity=30,Balanced=10,Gain=5)[profile]),ocs_backend=if(!ocs_enabled) NA_character_ else 'optiSel_adultcap32_v12',
   family_size_policy=if(!ocs_enabled) 'equal' else 'variable',
   target_offspring=if(early) 1500L else 500L,
   ocs_max_contribution=if(!ocs_enabled) NA_real_ else .10,
   speed_breeding=speed,scenario_id=stem,stringsAsFactors=FALSE)
 }
 rows <- list()
 for(first in c(9L,6L)) for(method in c('Pheno','Pedigree','GS')) rows[[length(rows)+1L]] <- make(method,first)
 for(method in c('Pedigree','GS')) rows[[length(rows)+1L]] <- make(method,4L)
 rows[[length(rows)+1L]] <- make('GS',2L)
 for(method in c('Pedigree','GS')) rows[[length(rows)+1L]] <- make(method,6L,TRUE)
 rows[[length(rows)+1L]] <- make('GS',2L,TRUE)
 for(first in c(6L,2L)) for(i in 1:3) rows[[length(rows)+1L]] <-
   make('GS',first,TRUE,ocs_enabled=TRUE,profile=c('Diversity','Balanced','Gain')[i])
 result <- do.call(rbind,rows)
 # Preserve existing folder identifiers for the 18 retained schemes.
 result$folder_name <- sprintf('%02d_%s',c(seq_len(12L),16:21),result$scenario_id)
 if(!is.null(root)) {
  for(i in seq_len(nrow(result))) {
   settings <- new.env(parent=baseenv())
   sys.source(file.path(root,result$folder_name[i],'GlobalParameters.R'),envir=settings)
   validate_ocs_settings(as.list(settings))
   for(field in c('selection_method','first_parental_use','parent_selection_stage',
                 'training_stage','early_selection','progeny_per_cross','ocs_enabled','ocs_target_Ne','speed_breeding')) {
    key <- if(field=='selection_method') 'method' else field
    if(!identical(settings[[field]],result[[key]][i])) stop('Scenario settings differ from the catalogue: ',result$folder_name[i],' / ',field)
   }
  }
  result$folder <- file.path(root,result$folder_name)
 }
 result
}

# Compatibility entry point used by existing launchers.
all_breeding_scenarios <- function(root=NULL) breeding_scenarios(root)

