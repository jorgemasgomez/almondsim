# Runtime monitoring

main.R runs all 25 schemes with four replicates and the bank in Haplotipos_2026-10-05_60000 (60,000 positions, 150 individuals).
Pheno/GS replicates use up to four parallel workers with six threads each.
Pedigree runs one replicate at a time in one licensed worker, after Pheno/GS.
The license is checked in a separate process before the experiment starts.

Results are outside Git, in almondsim_data/redesigned_runs.
almondsim_data/latest_run.txt identifies the current run.
Each rep_XXX contains progress_parallel.json / progress_pedigree.json,
execution logs, and timings_parallel.tsv / timings_pedigree.tsv.
Progress includes scenario, simulation year and update time.
Each scenario's results are saved immediately when that scenario completes.

Timings contain replicate, scenario, year, task, wall seconds, CPU seconds,
start timestamp and success. Tasks include founder initialization, phenotyping,
model fitting, predictions, selection, renewal, OCS/OPR, crossing and recording.
Times are inclusive: nested tasks overlap and must not be added together.
Use redesign_year rows for total yearly time and redesign_setup for initialization.
total_runtime.rds records complete experiment wall time, including I/O and coordination.

Memory-monitored configuration: up to four replicas at a time. Progress and timings include working-set, private and peak-private memory. diagnostics_*.log preserves warning messages, failing calls and the call stack. Timing-write errors do not abort a breeding simulation.
