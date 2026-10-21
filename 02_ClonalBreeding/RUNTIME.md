# Runtime and diagnostics

Progress: `experiment_progress.json`; each replicate has `progress_parallel.json`
and `progress_pedigree.json`, per-scheme status and errors, diagnostics and logs.
`timings_parallel.tsv` / `timings_pedigree.tsv` contain inclusive task times and memory.
Use `redesign_year` for annual duration and `redesign_setup` for initialization;
do not sum nested task times. `total_runtime.rds` records total elapsed wall time.

Allocation traces preserve validated native OCS plans, contributions, bounds and ages.
Result fields distinguish parentSelectionYear, firstUseOfSelectedParents and firstParentalUse.
The last year's queued selections are recorded but do not produce an extra future year.
