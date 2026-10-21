# Active manuscript breeding design

Years in scheme names denote **first parental use in crossing**.
The crossing/seed cycle is year 0; juvenile stages are 1â€“3; PT 4â€“5;
AT 6â€“8; ET 9â€“11. PT aggregate h2 is .20, AT .40, ET .50.
The fixed cuts are 500 -> 25 after PT and 25 -> 8 after AT.
Eight selected individuals completing AT in year 8 can first cross in year 9.
Parents completing PT in year 5 can first cross in year 6.
Predictive selection before own phenotyping permits first use in year 4.
Speed breeding permits pollen from year 2, with mothers from year 4.

Baseline replacement is 8 new + 24 retained. Increased selection intensity
uses 20 crosses x 75 offspring -> 500 selected seedlings in non-OCS arms.
All six OCS arms retain 1,500 offspring -> 500 selected seedlings, but allow variable
family sizes and any number of parental pairs and contributing trees. The ordinary
32-parent pool and fixed replacement rates continue to apply only to non-OCS arms.

OCS uses one optimizer variable per genotype and Ne-based controls 30/10/5
for Diversity/Balanced/Gain. GS6 retains the 10% individual cap; speed breeding
uses a 3.2% adult cap and 10% young cap, requiring at least 16 adult contributors. A single native
opticont max.BV fit uses the package-default solver. Native noffspring performs
integer rounding. Fixed role groups use matings; mixed-role speed breeding
uses a directed integer lpSolve allocation preserving total native counts.
There is no custom contribution repair or local refinement. Age, compatibility,
no selfing and cap remain hard constraints; rounded coancestry excess is warned
and audited. See [MANUAL_OCS.md](MANUAL_OCS.md) for the precise recurrence,
standing-parent reference and generation-interval proxy. Ne controls are not
guarantees of realized effective population size.

| Folder | First cross year | Field-evaluation years before first use | Early selection | OCS Ne control |
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
| 16_GS_Y6_Early1500_OCS_Diversity | 6 | 2 | TRUE | 30 |
| 17_GS_Y6_Early1500_OCS_Balanced | 6 | 2 | TRUE | 10 |
| 18_GS_Y6_Early1500_OCS_Gain | 6 | 2 | TRUE | 5 |
| 19_GS_SpeedBreeding_Y2_Early1500_OCS_Diversity | 2 | 0 | TRUE | 30 |
| 20_GS_SpeedBreeding_Y2_Early1500_OCS_Balanced | 2 | 0 | TRUE | 10 |
| 21_GS_SpeedBreeding_Y2_Early1500_OCS_Gain | 2 | 0 | TRUE | 5 |

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
Completed trial evaluations are generated with AlphaSimR::setPheno() at PT2,
AT3 and ET3, with reference h2 targets .20/.40/.50. These are separate trials,
evaluated in different calendar years, not overlapping cumulative means.
Intermediate stages retain the last completed evaluation; no yield phenotype
is available before PT2. The stage environment is P[year], without averaging
previous years. The same environmental calendar is shared within a replicate.
Each individual/stage record is cached and shared between trial and parent
copies. Retained parents keep their completed records, and acquire the next
stage evaluation according to age. No duplicate phenotype is generated.

Native residual variance equals A0*(1/h2-1)-GE0, with reps=1. With reference
A0=GE0=150000, PT/AT/elite residual variances are 450000/75000/0. Elite still
has GxE variation. Targets describe the reference population over random
environments, not necessarily realized heritability after selection.

GS-9 and Pedigree-9 use decision-specific repeated records. The PT cut uses
1500 original PT observations from the current PT cohort and its two predecessors.
AT selection and year-9 parental recruitment use all 1500 original PT observations
from the SAME three cohorts completing AT, plus the 75 end-of-AT observations of
their 25 survivors per cohort. These are 1575 observations on 1500 individuals;
75 IDs have paired PT/AT records with their own record years. Candidate own-record
membership, cohort matching and three-year record separation are asserted.

Both Y9 models retain fixed record Year and use the same independent-stage
relative non-genetic variance matrix diag(1, 0.375). Reference absolute
variances are PT 600000 and AT 225000, with zero cross-stage non-genetic
covariance. Genetic effects remain shared across observations of the same
individual. Overall residual scale and genetic/marker variance are estimated
separately in each fit. The former overlapping-error covariance is not used.
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
Parents retain completed-stage evaluations and receive the next evaluation according to age.
New candidates and existing parents are ranked separately for exact renewal.

Selections made after PT/AT are queued for next year's crosses. Before-phenotyping
arms recruit in the current year. Every crossing plan is validated against S-locus
and age restrictions. Non-OCS uses 20 equal-size crossing events; OCS produces
1,500 offspring distributed across variable-size families.
Repeated parental pairs are permitted. Rejected pollen does not change the parents
of an event or its allocated family size. Failures remain explicit.

## Version and provenance

The current full source design is `PT_AT_ET_stagePheno_v13`. The OCS optimization backend remains `optiSel_adultcap32_v12`.
New runs require a burn-in generated with the current source design.
Archived results retain their original frozen sources and versions. All 18 retained schemes belong to the manuscript comparison.
Early selection, cohort calendars and OCS constraints are unchanged. Stage phenotyping and the PT/AT error model differ from historical accumulated-evaluation runs; see [stage phenotyping](STAGE_PHENOTYPING.md).
For N=1500, a genotype contributes its combined maternal/paternal count divided
by 3000; the 10% cap allows at most 300 participations, and the 3.2% adult
cap in speed breeding allows at most 96. Audits preserve timings,
Ne control, parental reference, L proxy, counts, allocation mode and departures.
Updating local sources does not itself launch simulations or publish to GitHub.
