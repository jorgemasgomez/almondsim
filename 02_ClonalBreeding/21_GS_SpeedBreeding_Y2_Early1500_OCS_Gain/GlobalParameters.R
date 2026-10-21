selection_method <- "GS"
first_parental_use <- 2L
parent_phenotyping_years <- 0L
parent_selection_stage <- "Juvenile2"
parent_candidate_age <- 2L
training_stage <- "PT2"
selection_timing <- "before_own_phenotyping"
early_selection <- TRUE
progeny_per_cross <- 75L # defines total 20*75=1500; OCS family sizes are variable
seedlings_keep <- 500L
new_parents_per_year <- NA_integer_
renewal_percent <- NA_real_
nParents <- 32L # initial/burn-in pool only; not an OCS contributor limit
ocs_enabled <- TRUE
ocs_target_Ne <- 5
ocs_backend <- "optiSel_adultcap32_v12"
speed_breeding <- TRUE
minimum_mother_age <- 4L
minimum_pollen_age <- 2L

# OCS uses a fixed offspring total with variable family sizes and no parent-count limit.
ocs_max_contribution <- 0.10
