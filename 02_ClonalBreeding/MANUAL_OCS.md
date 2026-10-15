# Native optiSel offspring allocation (v9)

This file retains its historical name as a documentation entry point.

## Workflow

1. The current parents and eligible incoming candidates receive predictions from
   the same genomic model. Founder reference allele frequencies remain fixed.
2. Each eligible tree is represented by its maternal and/or paternal role.
   `optiSel::opticont` receives the role relationship matrix, breeding values,
   native individual bounds and native linear constraints limiting both copies
   of each tree jointly to 10%. Maternal and paternal totals are each one half.
3. Native `min.Kin` and `max.BV` fits define continuous endpoint coancestries.
   The target B=minK+alpha*(maxK-minK) uses alpha .10/.20/.40. A native `max.BV`
   fit uses this bound, adjusted only for optiSel's additive drift constant.
4. `optiSel::noffspring(..., N=1500, random=FALSE)` performs its own deterministic
   integer offspring rounding. Our code does not modify the resulting counts.
5. Rounded contributions and group coancestry are checked. A small coancestry
   excess is warned and recorded while native counts are preserved. The 10%
   combined-tree cap remains a hard check. No repair or retry is applied.
6. `optiSel::matings` pairs the fixed offspring margins. Its public solver hook
   passes the native objective and constraints to lpSolve and adds zero constraints
   for selfing, S-incompatible crosses and age-ineligible roles. There is no
   additional contribution objective or search in this adapter.
7. Each returned pair produces exactly its assigned number of offspring, with
   incompatible pollen rejected while its parents remain fixed.

There is no budget of 20 OCS pairs and no upper limit of 32 OCS trees. The number
of offspring remains 1,500 and early genomic selection retains 500. Non-OCS
scenarios retain their original family sizes, crossing budget and parent pool.

## Boundaries of the native workflow

Native contribution optimization does not prove that subsequent integer offspring
counts can always satisfy every biological pairing restriction. Native rounding
can also exceed the prescribed group-coancestry bound by a small amount. These
conditions are detected and reported, not concealed or privately repaired.
Coancestry rounding departures do not stop the scheme; infeasible biological
pairings and individual-cap violations do stop it.
The 10% limit combines both roles of the same tree; it is not 10% for each role.

Only biologically valid integer allocations satisfying the individual cap are
used to simulate offspring. A continuous coancestry target can be slightly exceeded
after native rounding; the excess is stored, not silently treated as compliance. With N offspring,
c_i=(maternal offspring_i+paternal offspring_i)/(2N). Pair coancestry is weighted
by each pair's offspring count. No singleton-dimension patch is needed in valid
v9 allocations: the cap requires at least five contributors in each role.

## Source and provenance

- `ocs_discrete.R`: genomic coancestry, compatibility and weighted metrics only.
- `ocs_optisel.R`: native optimization calls, restriction adapter and validation.
- `manual_ocs_adapter.R`: connects native allocations to the annual engine.
- `compatible_crosses.R`: creates the explicitly assigned variable-size families.
- `tests/check_ocs_native_offspring.R`: small allocation and integration fixtures;
  it does not launch breeding experiments.

Retired v8 helpers were preserved outside the repository. The design identifier is
`PT_AT_ET_optiSel_native_offspring_v9`. Existing frozen experiments are unchanged.
