# Optimal contribution selection

## Active settings

OCS is activated explicitly with `ocs_enabled <- TRUE`. The profile is defined
by `ocs_target_Ne`: **30 for Diversity, 10 for Balanced and 5 for Gain**.
Non-OCS schemes use `ocs_enabled <- FALSE` and `ocs_target_Ne <- NA_real_`.

These are target-Ne control parameters, not estimates or guarantees of realized
effective population size. Each decision uses the annual coancestry bound

    B = 1 - (1 - K_reference) * (1 - 1/(2*Ne))^(1/L).

`K_reference` is standing-parent mean genomic coancestry under equal genotype
weights. `L` is their mean age at the coming crossing decision, used as a
generation-interval proxy. Founder SNP reference frequencies remain fixed.

## Contributions and biological constraints

`optimize_contributions.R` constructs the candidate pool, breeding values, ages,
S-genotype groups and compatibility matrix. `ocs_optisel.R` optimizes one
contribution variable per genotype using a single native
`optiSel::opticont("max.BV", ..., solver="default")` call.

- GS6: combined maternal/paternal contribution <=10% per genotype.
- Speed breeding: contribution <=3.2% per adult (age >=4) and <=10% per younger
  candidate; pollen-only contributions jointly <=50%. At least 16 adults must
  therefore contribute, and the allocated offspring are checked accordingly.
- Mutually incompatible S-genotype groups: joint contribution <=50%.
  Sharing one S allele permits crossing through compatible pollen.
- Mothers: age >=4. Pollen donors: age >=2 in speed breeding, >=4 otherwise.
- No selfing and exactly 1,500 offspring, followed by early selection of 500.

Parent and pair counts are free. There is no fixed replacement rate or maximum
parent count for OCS. Family sizes may differ.

## Integer counts and matings

Native `optiSel::noffspring(random=FALSE)` supplies offspring counts.

1. With equal-role eligibility (GS6), maternal/paternal shares are rounded
   natively and `optiSel::matings` assigns compatible pairs.
2. If pollen-only candidates supply half of the total contribution, they form
   the paternal group and the remaining genotypes the maternal group. Native
   rounding supplies the fixed role counts before mating allocation.
3. Otherwise, native rounding supplies 3,000 total parental participations.
   A directed `lpSolve` formulation assigns adult roles and compatible pairs
   while preserving each genotype's exact count and the 1,500-offspring budget.

The mating adapters exclude forbidden edges and minimize offspring-weighted
pair coancestry. They require solver status 0 with a 30-second limit. Native
counts are not repaired, refined or changed to force feasibility.

## Validation and outputs

Continuous and allocated coancestry are checked against the original K matrix.
Relative excess up to 10% of the absolute bound plus 1e-7 is permitted and
audited; larger departures fail explicitly. This does not relax contribution
caps, ages, compatibility, no selfing, native solver validity or offspring totals.

Annual outputs record `ocsEnabled`, `ocsTargetNe`, the bound, coancestry, expected
merit and allocation departures. Allocation audits retain native fits, counts,
contributions, parental ages, SI-group limits, solver status and inclusive times.
Failed attempts remain recorded without automatic retry.

The executable source design is `PT_AT_ET_optiSel_adultcap32_v12`; new runs require
a matching burn-in. Existing frozen experiments and their outputs retain their
own source versions and are not modified when the repository is updated.
