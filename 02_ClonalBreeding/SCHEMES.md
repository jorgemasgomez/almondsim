# Active manuscript breeding design

Years in scheme names denote **first parental use in crossing**.
The crossing/seed cycle is year 0; juvenile stages are 1–3; PT 4–5;
AT 6–8; ET 9–11. PT aggregate h2 is .20, AT .40, ET .50.
The fixed cuts are 500 -> 25 after PT and 25 -> 8 after AT.
Eight selected individuals completing AT in year 8 can first cross in year 9.
Parents completing PT in year 5 can first cross in year 6.
Predictive selection before own phenotyping permits first use in year 4.
Speed breeding permits pollen from year 2, with mothers from year 4.

Baseline replacement is 8 new + 24 retained. GS–6 renewal rates 50/75/100%
use 500 individuals with **no early selection**. Increased selection intensity
uses 20 crosses x 75 offspring -> 500 selected seedlings in non-OCS arms.
All nine OCS arms retain 1,500 offspring -> 500 selected seedlings, but allow variable
family sizes and any number of parental pairs and contributing trees. The ordinary
32-parent pool and fixed replacement rates continue to apply only to non-OCS arms.

OCS v10 uses one optimizer variable per genotype, a combined 10% cap and
Ne-based controls 50/20/5 for Diversity/Balanced/Gain. A single native
opticont max.BV fit uses the package-default solver. Native noffspring performs
integer rounding. Fixed role groups use matings; mixed-role speed breeding
uses a directed integer lpSolve allocation preserving total native counts.
There is no custom contribution repair or local refinement. Age, compatibility,
no selfing and cap remain hard constraints; rounded coancestry excess is warned
and audited. See [MANUAL_OCS.md](MANUAL_OCS.md) for the precise recurrence,
standing-parent reference and generation-interval proxy. Ne controls are not
guarantees of realized effective population size.

| Folder | First cross year | Own phenotyping campaigns | Early selection | OCS Ne control |
|---|---:|---:|---|---:|
| 01_Pheno_Y9 | 9 | 5 | FALSE | NA |
| 02_Pedigree_Y9 | 9 | 5 | FALSE | NA |
| 03_GS_Y9 | 9 | 5 | FALSE | NA |
| 04_Pheno_Y6 | 6 | 2 | FALSE | NA |
| 05_Pedigree_Y6 | 6 | 2 | FALSE | NA |
| 06_GS_Y6 | 6 | 2 | FALSE | NA |
| 07_Pedigree_Y4 | 4 | 0 | FALSE | NA |
| 08_GS_Y4 | 4 | 0 | FALSE | NA |
| 09_GS_SpeedBreeding_Y2 | 2 | 0 | FALSE | NA |
| 10_Pedigree_Y6_Early1500 | 6 | 2 | TRUE | NA |
| 11_GS_Y6_Early1500 | 6 | 2 | TRUE | NA |
| 12_GS_SpeedBreeding_Y2_Early1500 | 2 | 0 | TRUE | NA |
| 13_GS_Y6_R50 | 6 | 2 | FALSE | NA |
| 14_GS_Y6_R75 | 6 | 2 | FALSE | NA |
| 15_GS_Y6_R100 | 6 | 2 | FALSE | NA |
| 16_GS_Y6_Early1500_OCS_Diversity | 6 | 2 | TRUE | 50 |
| 17_GS_Y6_Early1500_OCS_Balanced | 6 | 2 | TRUE | 20 |
| 18_GS_Y6_Early1500_OCS_Gain | 6 | 2 | TRUE | 5 |
| 19_GS_SpeedBreeding_Y2_Early1500_OCS_Diversity | 2 | 0 | TRUE | 50 |
| 20_GS_SpeedBreeding_Y2_Early1500_OCS_Balanced | 2 | 0 | TRUE | 20 |
| 21_GS_SpeedBreeding_Y2_Early1500_OCS_Gain | 2 | 0 | TRUE | 5 |
| 22_GS_Y4_Early1500_OCS_Diversity | 4 | 0 | TRUE | 50 |
| 23_GS_Y4_Early1500_OCS_Balanced | 4 | 0 | TRUE | 20 |
| 24_GS_Y4_Early1500_OCS_Gain | 4 | 0 | TRUE | 5 |

## Prediction and annual events

All GS and pedigree models include 1500 original end-of-PT records. Y9 AT
decisions additionally retain 75 paired AT records, as specified below.
Models are fitted separately for each decision, in EVERY predictive scheme:

* PT cut (500 -> 25): current complete PT2 cohort and the previous two PT2
  cohorts. The current candidates' own PT records MUST be included.
* AT cut (25 -> 8): full original PT2 populations of the current AT3 cohort
  and previous two AT3 cohorts. The AT survivors' own ORIGINAL PT records
  MUST be included. GS-9 and Pedigree-9 additionally use their 75 AT records;
  other schemes retain PT-only training.
* Parent recruitment: first-use year 9 uses the AT-anchored model; year 6
  uses the PT-anchored model. Before-phenotyping recruitment and early
  selection use the previous year's latest-PT model, because no own yield
  record exists yet.

Six full PT archives retain the three-year cohort lag. Missing/mismatched
archives or decision models without the candidates' own records fail explicitly.
Early1500 training uses the 500 retained seedlings per cohort, not all juveniles.
Parental retention uses the same model as recruitment; this update does not
add older parents' own records beyond the defined 1500-record window.

Annual environment probabilities are sampled as in the original program:
P[year] ~ Uniform(0,1), transformed by AlphaSimR's GxE trait to an annual normal
covariate. A separate seed per replicate makes the entire calendar common to
all paired scenarios, including burn-in, independently of scenario RNG draws.
Each evaluated individual's environmental covariates are averaged over its
phenotyped years (PT: years 4-5; end AT: years 4-8; end ET: years 4-11).
Aggregate residuals remain correlated; they are not stored as independent annual
phenotypes. Their marginal variance is calibrated against founder additive
variance: sigma_A0^2(1/h2-1) minus sigma_GxE0^2/nYears, avoiding double counting
annual GxE. This is a stage precision model, not raw repeated annual residuals.
The targets refer to the reference population, averaged over random annual
environments; realized precision may differ after selection and genetic change.
PT/AT/ET h2 assumptions remain .20/.40/.50. Existing parents continue accumulating
environmental evaluations; no duplicate evaluation is counted within a year.
No additional, independent common Year offset is invented: the original GxE
covariate itself changes yearly genotype responses and can shift population means.

GS-9 and Pedigree-9 use decision-specific repeated records. The PT cut uses
1500 original PT observations from the current PT cohort and its two predecessors.
AT selection and year-9 parental recruitment use all 1500 original PT observations
from the SAME three cohorts completing AT, plus the 75 end-of-AT observations of
their 25 survivors per cohort. These are 1575 observations on 1500 individuals;
75 IDs have paired PT/AT records with their own record years. Candidate own-record
membership, cohort matching and three-year record separation are asserted.

Both Y9 models retain fixed record Year and use the same design-calibrated relative
residual matrix [[1, 0.3416666667], [0.3416666667, 0.375]] (rho=0.5579393),
including the aggregate-error recursion and expected shared GxE. Reference absolute
values are PT 600000, AT 225000, covariance 205000. Overall residual scale and
genetic/marker variance are estimated separately in each fit. This is a marginal
working covariance calibrated from the simulated design, not an estimate based
on true individual genetic values or the complete realised GxE covariance.
GS-9 uses native AlphaSimR::solveRRBLUP after whitening response, fixed design and
markers; Pedigree-9 uses ASReml with proportional residual constraints. Its PT-only
fit keeps the common residual. Both stages preserve fixed Year.

Other GS schemes keep AlphaSimR::RRBLUP and 1500 PT-only training observations;
other pedigree schemes keep ASReml with a common PT residual. Decision-specific
cohorts and own-record inclusion remain mandatory at PT and AT cuts. No new AT
records are added to years 6, 4 or 2. No additional parental records outside the
three training cohorts are appended. OCS uses the validated optiSel contribution and mating backend described above.

The PT-anchored model updated after the year's evaluations is carried
to the next year's early selection and juvenile parental selection.
Parents continue to be evaluated according to age. Aggregate precision is calibrated jointly with the annual GxE contribution.
New candidates and existing parents are ranked separately for exact renewal.

Selections made after PT/AT are queued for next year's crosses. Before-phenotyping
arms recruit in the current year. Every crossing plan is validated against S-locus
and age restrictions. Non-OCS uses 20 equal-size crossing events; OCS produces
1,500 offspring distributed across variable-size families.
Repeated parental pairs are permitted. Rejected pollen does not change the parents
of an event or its allocated family size. Failures remain explicit.

## Version and provenance

The current OCS design is `PT_AT_ET_optiSel_single_genotype_Ne_v10`.
Earlier matched burn-ins can be reused with their seeds and environments;
historical future results must retain their original versions. Only GS6 and
speed-GS2 OCS belong to the paper; exploratory GS4 OCS is not launched.
The early selection, calendars, phenotyping and prediction models are unchanged.
For N=1500, a genotype contributes its combined maternal/paternal count divided
by 3000; the 10% cap allows at most 300 participations. Audits preserve timings,
Ne control, parental reference, L proxy, counts, allocation mode and departures.
Updating local sources does not itself launch simulations or publish to GitHub.
