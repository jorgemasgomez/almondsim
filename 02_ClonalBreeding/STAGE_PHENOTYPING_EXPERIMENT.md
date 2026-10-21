# Experimental stage phenotyping

This branch changes the scientific phenotype model. It must not be merged into
the main pipeline merely because a single-replicate pilot completes.

## Separate PT, AT and elite trials

`update_phenotypes()` calls `AlphaSimR::setPheno()` at PT2, AT3 and ET3 only.
PT and AT represent distinct trials evaluated in different calendar years.
Intermediate years retain the most recently completed trial record; before PT2
there is no available yield phenotype. Parent and trial copies of a genotype
share a cached record, preventing duplicate evaluation with a different error.
Parents retain their completed evaluation rather than being re-phenotyped
independently every year. A genotype recruited before trial selection may receive
its next completed-stage evaluation as a parent even if its trial copy was cut.

Each call uses the pre-generated calendar environment `p = environment_p[year]`.
There is no average of previous environments. Native AlphaSimR GxE sensitivities
remain fixed for each genotype; residual draws of distinct stage trials are
independent. Evaluations still share the genotype's genetic value and GxE slope.

## Variance calibration

The aggregate target non-genetic variance is `A * (1 / h2 - 1)`, with the original
reference genetic variance `A`. A single stage environment contributes the
reference GxE variance `GE`, so native `setPheno()` receives
`varE = A * (1 / h2 - 1) - GE` and `reps = 1`. These are entry-mean errors, not
plot-level errors or explicit annual observations.

With A = GE = 150,000, residual budgets are 450,000 (PT, h2 = .20),
75,000 (AT, h2 = .40) and zero (elite, h2 = .50). Elite still has GxE variation;
zero here concerns the additional residual only. Higher target heritability
would be infeasible without changing the assumed stage GxE variance.

Year-9 genomic and pedigree models retain 1,500 original PT observations and
75 AT observations for the same cohort identities. The non-genetic PT/AT
working covariance is now zero, and the relative variance matrix is
`diag(c(1, .375))`. Genomic whitening and the ASReml `idh(Stage)` model use the
same shape. The legacy accumulated-error covariance must not be reused.

## Pilot and validation

The pilot requests one replicate of each of the 18 current schemes, 40 future
years each, from one NEW shared 40-year burn-in after pipeline filling.
Historical burn-ins are incompatible. The founder bank, founder sampling seed,
calendar seed, trial sizes, parental rules and OCS settings remain unchanged.
Different random-number consumption and selection during burn-in mean that
this is a same-seed comparison, not an identical starting-population comparison.

Real AlphaSimR tests compare native `setPheno()` values exactly, check record
caching and archives, and verify independent PT/AT covariance. A real ASReml
preflight and 1,500 PT + 75 AT genomic/pedigree fitting preflights also passed.

Renaming-only regression tests passed on main before this scientific change.
Their historic phenotype/calibration signatures intentionally do not describe
this experimental branch. `tests/check_stage_phenotyping.R` tests the new model.

The pilot must report coverage, failures, annual accuracy, gain relative to its
own burn-in, final genetic variance, and wall time. A single replicate assesses
function and the magnitude of changes; it cannot establish equivalence or
preservation of the scheme ranking. No automatic failed-scheme retries.
