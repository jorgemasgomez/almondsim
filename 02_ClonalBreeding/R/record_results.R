# Record annual outcomes by population stage and combine records with compatible columns.
# Function names describe each operation; all launchers use this shared implementation.

# Record population sizes, genetic means, genetic variances and selection accuracies for the current year.
summarize_annual_populations <- function(state,scenario,rep,accuracy,entrants,parent_stage) {
 row <- list(year=state$year,rep=rep,scenario=scenario,accSel=accuracy,
             accParents=state$parent_accuracy,newParents=entrants,parentSelectionStage=parent_stage)
 for(stage in c('Parents',state$config$stages)) {
  pop <- if(stage=='Parents') state$parents else state$pop[[stage]]
  row[[paste0('n_',stage)]] <- if(is.null(pop)) 0L else nInd(pop)
  row[[paste0('meanG_',stage)]] <- if(is.null(pop)) NA_real_ else as.numeric(meanG(pop))[1L]
  row[[paste0('varG_',stage)]] <- if(is.null(pop)) NA_real_ else as.numeric(varG(pop))[1L]
 }
 row$meanG <- row$meanG_Seedlings;row$varG <- row$varG_Seedlings
 haplo <- pullMarkerHaplo(state$pop$Seedlings,markers=locuscompt_position_vector,haplo='all',simParam=SP)
 row$allelesSI <- length(unique(apply(haplo,1,paste,collapse='_')))
 as.data.frame(row,stringsAsFactors=FALSE)
}

# Combine annual/status frames while preserving all columns when some records lack a field.
bind_result_rows <- function(frames) {
 frames<-Filter(function(x) !is.null(x)&&nrow(x)>0L,frames)
 if(!length(frames)) return(NULL)
 columns<-unique(unlist(lapply(frames,names)))
 do.call(rbind,lapply(frames,function(x) {
   for(column in setdiff(columns,names(x))) x[[column]]<-NA
   x[columns]
 }))
}


# Technical output fields are kept here so advance_year.R follows the biology.
# These inputs are the decisions already made by the annual cycle. Recording
# does not choose parents, fit a new model or perform additional selection.
record_annual_results <- function(state, config, method, first, rep, accuracy,
                                entrants, target, allocation, next_selection,
                                ocs_solution, candidates, seedlings_before) {
 scenario <- if(!is.null(config$scenario_id)) config$scenario_id else paste0(method,'_Y',first)
 row <- summarize_annual_populations(state,scenario,rep,accuracy,entrants,target)
 row$ocsMaxContribution <- if(!isTRUE(config$ocs_enabled)) NA_real_ else config$ocs_max_contribution
 row$firstParentalUse <- first
 row$parentCandidateAge <- if(first>=6L) first-1L else first
 row$parentPhenotypingYears <- if(first==9L) 5L else if(first==6L) 2L else 0L
 row$parentSelectionYear <- if(is.null(allocation)) NA_integer_ else allocation$selected_year
 row$firstUseOfSelectedParents <- if(is.null(allocation)) NA_integer_ else allocation$use_year
 row$selectedForNextYear <- if(is.null(next_selection)) 0L else next_selection$new
 row$trialCutStage <- 'PT2';row$eliteCutStage <- 'AT3'
 row$predictionSolver <- if(method=='GS') if(first==9L) 'AlphaSimR::solveRRBLUP' else 'AlphaSimR::RRBLUP' else if(method=='Pedigree') 'ASReml' else NA_character_
 row$trainingPTRecords <- if(method=='Pheno') NA_integer_ else sum(build_training_population(state,candidates)@misc$training_record_stage=='PT')
 row$trainingATRecords <- if(method=='Pheno') NA_integer_ else if(first==9L) 75L else 0L
 row$trainingCohortAnchor <- if(method=='Pheno') NA_character_ else state$training_cohort_stage
 row$annualEnvironmentP <- state$environment_p[[as.character(state$year)]]
 row$trainingStage <- if(method=='Pheno') NA_character_ else state$training_stage
 row$trialTrainingAnchor <- if(method=='Pheno') NA_character_ else 'PT2'
 row$eliteTrainingAnchor <- if(method=='Pheno') NA_character_ else 'AT3'
 row$trialTrainingPTRecords <- if(method=='Pheno') NA_integer_ else 1500L
 row$eliteTrainingPTRecords <- if(method=='Pheno') NA_integer_ else 1500L
 row$trialTrainingATRecords <- if(method=='Pheno') NA_integer_ else 0L
 row$eliteTrainingATRecords <- if(method=='Pheno') NA_integer_ else if(first==9L) 75L else 0L
 row$trainingRecords <- if(method=='Pheno') NA_integer_ else if(first==9L) 1575L else 1500L
 row$trainingVariant <- if(method=='Pheno') NA_character_ else if(first==9L) 'Y9_PT1500_plus_AT75_calibrated' else 'decision_specific_PT1500'
 row$residualShape <- if(method=='Pheno') NA_character_ else if(first==9L) 'design_calibrated_PT_AT' else 'common_PT'
 row$ownPTRecordsValidated <- method!='Pheno'
 row$trainingCohorts <- if(method=='Pheno') NA_integer_ else min(state$config$training_cohorts,length(state$stage_history[[state$training_cohort_stage]]))
 row$trainingIndividuals <- if(method=='Pheno') NA_integer_ else nInd(build_training_population(state,candidates))
 row$earlySelection <- isTRUE(config$early_selection)
 row$progenyPerCross <- if(is.null(ocs_solution)) config$progeny else NA_integer_
 row$variableFamilySizes <- !is.null(ocs_solution)
 row$numberOfParentalPairs <- if(is.null(state$speed_cross_trace)) NA_integer_ else state$speed_cross_trace$distinct_pair_count
 row$maxContributingParents <- if(!isTRUE(config$ocs_enabled)) config$parents else NA_integer_
 row$seedlingsBeforeSelection <- seedlings_before
 row$seedlingsAfterSelection <- nInd(state$pop$Seedlings)
 row$targetNewParents <- if(!isTRUE(config$ocs_enabled)) config$entrants else NA_integer_
 row$targetRenewalPercent <- if(!isTRUE(config$ocs_enabled)) 100*config$entrants/config$parents else NA_real_
 row$realizedRenewalPercent <- 100*entrants/nInd(state$parents)
 row$renewalRule <- if(!isTRUE(config$ocs_enabled)) 'exact' else 'OCS_free'
 row$parentsUsed <- nInd(state$parents)
 row$meanParentalAge <- mean(state$year-state$birth[as.character(state$parents@id)])
 trace <- state$speed_cross_trace
 row$minimumPollenAge <- config$minimum_pollen_age;row$minimumSeedParentAge <- config$minimum_mother_age
 row$meanMotherAge <- mean(trace$mother_ages);row$meanFatherAge <- mean(trace$father_ages)
 row$weightedParentalAge <- (row$meanMotherAge+row$meanFatherAge)/2
 row$pollenFromAge2 <- mean(trace$father_ages==2L)
 row$pollenFromUnder4 <- mean(trace$father_ages<4L)
 row$ocsEnabled <- isTRUE(config$ocs_enabled);row$ocsTargetNe <- config$ocs_target_Ne;row$ocsImplementation <- if(!isTRUE(config$ocs_enabled)) NA_character_ else 'optiSel_adultcap32_v12'
 row$ocsBound <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$bound
 row$ocsGenomicCoancestry <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$metrics$kinship
 row$ocsMateCoancestry <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$metrics$mate_kinship
 row$ocsExpectedBV <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$metrics$gain
 row$ocsEffectiveParents <- if(is.null(ocs_solution)) NA_real_ else 1/sum(ocs_solution$metrics$contribution^2)
 row$ocsBoundSatisfied <- if(is.null(ocs_solution)) NA else ocs_solution$metrics$kinship<=ocs_solution$bound+1e-8
 row$ocsCoancestryExcess <- if(is.null(ocs_solution)) NA_real_ else max(0,ocs_solution$metrics$kinship-ocs_solution$bound)
 row$ocsCoancestryExcessRelative <- if(is.null(ocs_solution)) NA_real_ else ocs_solution$solver_details$optimization_audit$rounded_kinship_excess_relative
 row$designVersion <- config$design_version
 row
}
