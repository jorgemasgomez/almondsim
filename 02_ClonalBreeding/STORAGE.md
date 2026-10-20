# Code and external data

The active catalogue contains the 18 manuscript schemes. Run `main.R` explicitly for all schemes,
or `<numbered scenario>/00RUNME.R` for one. Sourcing the engine only defines functions.
Default: 10 replicates, 40 burn-in years and 40 future years, preceded by 14 fixed-founder
initialization cycles. Full successful future coverage is 7,200 unique rows (18 x 10 x 40).

New runs are written outside this repository to the sibling
`almondsim_data/manuscript_runs` directory. `ALMONDSIM_DATA_DIR` changes the external
data directory; `ALMONDSIM_FOUNDER_BANK` changes the bank150/60000 input directory.
The data directory must remain outside the source repository. All snapshots, haplotypes,
logs, results and plots remain external. A run freezes its source in its own `code` folder.

Pheno/GS use four workers with six threads; pedigree runs sequentially in one worker,
ASReml workspace 16gb. A failed scenario stops only that scenario in that replicate.
Other scenarios and replicates continue. Plot coverage explicitly excludes incomplete
trajectories. Historical code is backed up outside this source tree during migration.

Current source design: PT_AT_ET_optiSel_adultcap32_v12. OCS uses native offspring
allocation with variable family sizes, no upper limit on parent/pair count, and a
10% combined contribution cap in GS6. Speed breeding uses a 3.2% adult cap
and a 10% young cap, requiring at least 16 contributing adults. Frozen earlier experiments are not modified.
