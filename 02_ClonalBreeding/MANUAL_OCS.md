# Single-genotype optimal contributions (v10)

## Selection and constraints

The continuous optimizer has one variable per genotype, not separate maternal
and paternal copies. One `optiSel::opticont("max.BV", ..., solver="default")`
call maximizes predicted genetic merit. The combined individual cap is 10%.
Candidates that can supply pollen only have a joint upper bound of 50%; an
analogous bound applies to any mother-only group. These aggregate bounds are
necessary but do not guarantee compatible integer matings.

Diversity, Balanced and Gain use control parameters Ne = 50, 20 and 5,
respectively. The annual genomic coancestry bound is

    B = 1 - (1 - K_reference) * (1 - 1/(2*Ne))^(1/L).

K_reference is the mean genomic coancestry of the standing parental pool,
with equal genotype weights. L is the mean age of that pool at the coming
crossing decision, used as an explicit generation-interval proxy. These are
Ne-based control parameters, not estimates or guarantees of realized effective
population size. The matrix is founder-centred genomic coancestry; the target
must be interpreted on that reference scale. Audit files preserve both inputs.
There are no minimum-coancestry/maximum-merit endpoint fits. Legacy `ocs_alpha`
fields remain internal scheme dispatch markers only; alpha no longer defines
the bound. Annual outputs record ocsAlpha=NA and ocsTargetNe explicitly.

## Native rounding and biological allocation

`optiSel::noffspring(random=FALSE)` supplies integer counts without custom
repair or local contribution refinement:

* When all candidates are eligible in both roles (GS6), equal continuous
  maternal/paternal shares are rounded natively; `optiSel::matings` allocates
  fixed role counts with a forbidden-edge lpSolve adapter.
* When pollen-only candidates saturate half of the contribution (within 1e-7),
  they form the paternal group and the remaining genotypes the maternal group.
  Native group rounding preserves 1500 counts in each role before matings.
* Otherwise native noffspring allocates 3000 total parental participations.
  A single bookkeeping group uses oc=c/2 and N=3000; no biological sex is
  assigned by that bookkeeping label. A directed integer lpSolve model jointly
  assigns roles and pairs. Each genotype's total participation count is fixed,
  each forbidden biological edge is excluded, and 1500 offspring are required.
  The objective minimizes offspring-weighted parental-pair coancestry.
  This general-case allocator is custom glue and is NOT optiSel::matings.

All allocations enforce maternal age >=4, paternal age >=2 in speed breeding
(otherwise >=4), gametophytic compatibility, no selfing, the offspring budget
and the combined 10% cap. Families can differ in size; parent/pair count is free.
No native integer counts are changed to force feasibility. Infeasible counts,
cap violations or unsuccessful solvers fail that scenario/replicate without
retry. Pair allocation has a 30-second limit and requires lpSolve status 0.

Coancestry is checked against the original matrix. The user-authorized policy
of18October accepts continuous and rounded relative excesses up to10% (with
1e-7 absolute numerical tolerance), recording warnings and absolute/relative
audits. Larger departures fail explicitly. Native solver validity remains required. Rounded excess does not relax age, compatibility or
cap checks. No claim of a globally optimal discrete OCS solution is made.

## Provenance

Design: `PT_AT_ET_optiSel_single_genotype_Ne_v10`. Old v7/v8/v9 future results
retain their original code and are excluded from v10 OCS comparisons. Saved
matched burn-ins, seeds and environments can be reused without relabelling
historical future trajectories. Non-OCS/Pedigree models and crossing rules are
unchanged. The paper includes GS6 and speed-GS2 OCS only, not exploratory GS4.
Audits include optimizer size, solver, continuous and realized contributions,
allocation mode, bounds, Ne target, reference coancestry, L proxy and timings.
