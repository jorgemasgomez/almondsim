# OCS backend (compatibility filename)

`ocs_discrete.R` provides genomic coancestry, biological compatibility, plan metrics
and capped local refinement; `ocs_cap10_support.R` constructs the capped endpoints.
`manual_ocs_adapter.R` combines current parents and eligible new
candidates, predicts all with the same genomic model, and applies directed S-locus
compatibility and maternal/paternal age masks before optimization.

Each parental appearance contributes 1/(2*20). Plans have 20 crossing events and at
most 32 contributing parents. The bound is minK + alpha*(maxK-minK) using feasible
capped endpoints. Diversity/Balanced/Gain use alpha .10/.20/.40. Contributions are
optimized with optiSel::opticont, rounded and refined under hard constraints;
optiSel::matings then minimizes pair coancestry conditional on role counts.
Optimization randomness does not consume the breeding stream. There is no automatic
replacement of failed scenarios or substitution of an uncapped endpoint plan.

The search is heuristic. A failed scheme is reported and isolated; all other schemes
continue. Endpoint collapse is recorded for interpretation. New parental slots are
not capped at eight in OCS; that quota applies only to exact-renewal comparisons.


## OCS contribution cap (v8, prepared for subsequent runs)

All OCS profiles now enforce a maximum combined maternal+paternal genetic contribution of 0.10 per tree. With 20 equal-size crosses this permits at most four of 40 parental participations, hence at least ten distinct contributors overall and at least five in each role. Maximum32 parents, age eligibility, S compatibility and no-selfing remain unchanged. The cap is enforced in continuous contributions, integer recovery and local refinement, with independent final validation.

Both alpha reference endpoints are recomputed under this cap; their values must not be mixed with the uncapped v7 endpoints. Integer endpoint construction uses a bounded shortlist and local refinement and is not claimed globally optimal. An infeasible capped design fails explicitly rather than substituting an uncapped plan.

A local copy of optiSel::matings retains matrix dimensions (`drop=FALSE`) for singleton maternal or paternal roles. This fixes the reproduced v7 failure without changing the installed package.

No simulation is launched by this source update. The running experiment uses its own frozen v7 code; its results have no 10% cap and remain labelled v7. New runs default to PT_AT_ET_optiSel_cap10_v8 and record ocsMaxContribution.
