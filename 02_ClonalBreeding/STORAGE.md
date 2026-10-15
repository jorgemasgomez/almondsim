# Code and external data

The active catalogue has 24 schemes; the manuscript excludes the three exploratory GS4 OCS schemes (21 schemes). Run `main.R` explicitly for all schemes,
or `<numbered scenario>/00RUNME.R` for one. Sourcing the engine only defines functions.
Default: 10 replicates, 40 burn-in years and 40 future years, preceded by 14 fixed-founder
initialization cycles. Full successful future coverage is 9,600 unique rows for the full catalogue (24 x 10 x 40), or 8,400 rows for the 21 manuscript schemes.

New runs are written outside this repository to the sibling
`almondsim_data/manuscript_runs` directory. `ALMONDSIM_DATA_DIR` changes the external
data directory; `ALMONDSIM_FOUNDER_BANK` changes the bank150/60000 input directory.
The data directory must remain outside the source repository. All snapshots, haplotypes,
logs, results and plots remain external. A run freezes its source in its own `code` folder.

Pheno/GS use four workers with six threads; pedigree runs sequentially in one worker,
ASReml workspace 16gb. A failed scenario stops only that scenario in that replicate.
Other scenarios and replicates continue. Plot coverage explicitly excludes incomplete
trajectories. Historical code is backed up outside this source tree during migration.
