# Code and external simulation data

The active experiment contains 25 scenarios, including GS-2 early1500 OPR Diversity (5%), Balanced (10%), and Gain (20%). See `02_ClonalBreeding/OPR_design.md`.

Run `02_ClonalBreeding/main.R` for the full comparison, or a scenario's `00RUNME.R` directly under `02_ClonalBreeding/<scenario>`. Results and burn-in snapshots for new active runs are generated outside this repository in the sibling `almondsim_data/redesigned_runs` directory. Set the environment variable `ALMONDSIM_DATA_DIR` to choose a different external directory. Legacy haplotypes are read from `<data directory>/Haplotypes`.

The 150-individual founder bank is an external input; set `config$founder_bank` to its directory. Generation is separate from this repository and the loader requires a complete manifest. No haplotypes, plots, results, burn-in snapshots, logs, or simulation run directories belong in Git.

Historical numbered scenario folders are archived outside the repository in `almondsim_data/archived_code`. The common `00_Burn_in` initialization remains a runtime dependency. Use the active launcher for the redesigned experiment. A pilot that was already running when storage changed retains its original directory until completion; it is excluded from Git.

Existing binary files removed from the current Git tree remain recoverable in earlier commits. This update does not rewrite repository history.
