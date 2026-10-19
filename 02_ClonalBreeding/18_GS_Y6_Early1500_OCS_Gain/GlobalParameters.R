selection_method <- "GS"
first_parental_use <- 6L
parent_phenotyping_years <- 2L
parent_selection_stage <- "PT2"
parent_candidate_age <- 5L
training_stage <- "PT2"
selection_timing <- "end_of_previous_year"
early_selection <- TRUE
progeny_per_cross <- 75L # defines total 20*75=1500; OCS family sizes are variable
seedlings_keep <- 500L
new_parents_per_year <- NA_integer_
renewal_percent <- NA_real_
nParents <- 32L # initial/burn-in pool only; not an OCS contributor limit
# Legacy dispatch marker only; no alpha interpolation in v10.
ocs_alpha <- 0.4
ocs_target_Ne <- 5
ocs_backend <- "optiSel_single_genotype_Ne_v10"
speed_breeding <- FALSE
minimum_mother_age <- 4L
minimum_pollen_age <- 4L

# OCS uses a fixed offspring total with variable family sizes and no parent-count limit.
ocs_family_sizes <- "variable"
ocs_max_contribution <- 0.10
ocs_max_parents <- NA_integer_
