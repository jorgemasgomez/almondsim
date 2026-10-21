# Parental renewal throughout one breeding year

`state$parents` holds the parents active for crossing this year. `state$pending`
holds a complete allocation chosen last year for activation this year. Candidate
populations in `state$pop` are separate cohorts, not the active parental pool.

## First parental use at age 6 or 9

1. At the beginning of the year, `advance_breeding_year()` checks that a pending
   allocation has `use_year == state$year` and activates its parents.
2. Those active parents produce this year's F1 offspring.
3. Trial populations are evaluated and their complete records archived.
4. `fit_decision_model()` prepares the PT model for age-6 use or the PT/AT
   model for age-9 use.
5. `select_parent_allocation()` selects the next parental pool: candidates are
   aged 5 in PT2 or aged 8 in AT3. The allocation is stored in `state$pending`.
6. Next year, at age 6 or 9, the selected individuals become active parents.

At the beginning of the first future year, the common burn-in parents remain
active and the burn-in's pending allocation is cleared. Each new scheme first
chooses its own pending allocation at that year's end.

## First parental use at age 2 or 4

After cohort advancement, the age-2 Juvenile2 or age-4 PT1 candidates are scored
with the previously available PT model, before their current-year phenotyping.
`select_parent_allocation()` replaces the active parents immediately, and the
new allocation is crossed in the same year. There is no one-year pending delay.
Speed breeding permits young candidates to supply pollen only; mothers must be
at least four years old.

## Ordinary replacement versus OCS

Without OCS, `renew_parent_pool()` calls `select_best_individuals()`, which uses
AlphaSimR's `selectInd()`. It recruits the best eight distinct incoming
candidates and retains the best 24 existing parents. These are two independent
rankings, not selection of the best 32 from a combined pool. Recruitment of
eight newcomers is obligatory even if they score below existing parents.

When parent recruitment and trial selection use identical scores, the eight
incoming PT parents are among the 25 individuals retained for AT. At AT3, the
eight new parents coincide with the eight elite selections under the same
ranking. Crossing remains a separate operation from selecting these parents.

With OCS, `optimize_parent_contributions()` instead optimizes over existing
parents and eligible incoming candidates. The number of contributors is free,
subject to the configured constraints. All candidates must remain available
before trial discarding because a useful OCS contributor need not be among
the highest-ranking trial selections by breeding value alone.

## Code map

| Operation | Implementation |
| --- | --- |
| Activate pending parents; cross; queue next year's allocation | `R/advance_year.R` |
| Score candidates; dispatch ordinary replacement or OCS; label use year | `R/select_parents.R` |
| Select the best individuals with AlphaSimR | `R/apply_selection.R` |
| Optimize OCS contributions and parental roles | `optimize_contributions.R` and `ocs_optisel.R` |

Renaming these functions does not change selection, populations, scores or RNG
use. A separate experimental branch changes phenotyping; that is a scientific
model change rather than part of this naming refactor.
