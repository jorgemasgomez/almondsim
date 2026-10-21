# Reviewing the breeding simulation

Start with [PIPELINE.md](PIPELINE.md) for the script connections and annual order.
All paths below are relative to `02_ClonalBreeding/`.

## Suggested reading order

1. **Choose a scheme:** `R/simulation_parameters.R` and `R/run_interface.R`.
   The former contains the scientific catalogue; the latter maps manuscript
   labels such as `GS6-E` onto stable scenario IDs. `GlobalParameters.R` in each
   scheme directory remains a catalogue validation input.
2. **Follow one year:** `R/advance_year.R`, starting at numbered comments 1-9.
   Check whether each decision uses previous or current evaluations and whether
   newly chosen parents are used immediately or in the following year.
3. **Check phenotypes:** `R/simulate_phenotypes.R`. Review annual environmental
   sharing, accumulated evaluation years and correlated aggregate errors.
4. **Check prediction:** `R/build_training_population.R`, then
   `R/fit_prediction_model.R` and `R/predict_genetic_values.R`. Check cohort
   identities and candidate own records, not only training-set sizes.
5. **Check selection:** `R/apply_selection.R` and `R/select_parents.R`.
6. **Check crosses/OCS:** `compatible_crosses.R`, `ocs_discrete.R`,
   `optimize_contributions.R` and `ocs_optisel.R`. See the dedicated OCS diagram.
7. **Check measurement and replication:** `R/record_results.R`,
   `R/run_replicate.R` and `R/run_experiment.R`.

## Questions and their code locations

| Review question | Read here |
| --- | --- |
| What differs between PH9, GS6-E and SB2-OG? | `R/simulation_parameters.R`: `redesign_scenarios()` |
| When are cohorts aged and when are parents crossed? | `R/advance_year.R`, `R/advance_cohorts.R` |
| Is seedling accuracy measured before early selection? | `R/advance_year.R`: step 3 |
| Which observations train each selection decision? | `R/build_training_population.R` |
| Where is the calibrated repeated PT/AT covariance applied? | `R/fit_prediction_model.R`: `redesign_calibrated_residual()`, `redesign_native_repeated_gs()`, `redesign_asreml_calibrated()` |
| Are PT and AT cuts made before or after parental decisions? | `R/advance_year.R`: steps 8-9; `R/apply_selection.R` |
| How is the target-Ne bound constructed and checked? | `optimize_contributions.R`, `ocs_optisel.R` |
| How are S alleles, parental ages and contribution caps enforced? | `optimize_contributions.R`, `ocs_discrete.R`, `ocs_optisel.R`, `compatible_crosses.R` |
| How are shared burn-ins, RNG streams and environmental histories maintained? | `R/run_experiment.R`, `R/run_replicate.R`, `R/simulate_phenotypes.R` |
| What happens after a scenario fails? | `R/run_replicate.R`: outcome/status handling |
| Is a running experiment protected from later source edits? | `R/freeze_pipeline_source.R` and the run's `code/source_manifest.csv` |
| Which population defines annual gain and variance? | `R/record_results.R`; derived burn-in-relative gain is calculated during post-processing |

## Running and inspecting a scheme

Run commands from the repository root with dependencies installed and the
external founder/data paths configured as described in the README.

```sh
# Information only: neither command starts a simulation.
Rscript 02_ClonalBreeding/run_scheme.R --help
Rscript 02_ClonalBreeding/run_scheme.R --list

# Explicit new experiment; never run this to inspect an old result.
Rscript 02_ClonalBreeding/run_scheme.R --scheme=GS6-E --replicates=1 --workers=1 --threads=6

# Selected schemes, with the same engine and annual sequence.
Rscript 02_ClonalBreeding/run_experiment.R --schemes=GS6-E,SB2-E --replicates=10 --workers=4 --threads=6
```

For an interactive R session:

```r
root <- normalizePath('02_ClonalBreeding', winslash='/')
source(file.path(root, 'redesigned_schemes.R')) # definitions only
config <- redesign_config()
config$reps <- 1L
config$workers <- 1L
config$threads <- 6L

select_breeding_schemes(root, 'GS6-E') # inspect settings; no run
# Call only when a new simulation is intended:
# result <- run_breeding_scheme(root, 'GS6-E', config)
```

## Refactor validation

```sh
Rscript 02_ClonalBreeding/tests/check_readable_pipeline.R
Rscript 02_ClonalBreeding/tests/check_ne_configuration.R
Rscript 02_ClonalBreeding/tests/check_ocs_native_offspring.R
```

The first check compares the executable syntax of 39 unchanged functions
against signatures captured from commit
`da0f9c8c2e7b46e10b972f84e4690229549d882e` before reorganization. It also checks
unique module definitions, scheme aliases, cohort movement, CLI dispatch,
source-freeze checksums and module loading in a separate R worker.
Comments, whitespace and function source locations do not affect the signatures.
The experiment controller has one intentional operational change: freezing the
new `R/` directory in addition to the existing top-level files. Its worker
scheduling and scientific settings remain unchanged. A separate expected
controller signature checks that only this source-copy block was replaced.

The Ne configuration check covers the six intentional configuration/dispatch
changes separately: an explicit OCS switch, target-Ne settings, removal of
unused output fields, and loading only the current source design. Catalogue
settings are checked against every scheme directory. `R/simulation_parameters.R`
is the executable catalogue; there is no separately maintained JSON copy.

The optional `analysis/export_pedigree.R` utility exports a saved pedigree to
CSV. It does not participate in simulation execution.

The final check exercises the current OCS allocation and crossing fixtures.
These checks do not launch breeding replicates and do not require an ASReml
licence. They do not claim to be a new validation of full 40-year trajectories
or of licensed ASReml numerical fits. The sources and results used for the
manuscript remain in their existing frozen runs, outside this refactor.
