# Almond breeding simulations with AlphaSimR

## Introduction

This repository contains AlmondSim scripts for almond breeding-program simulations using AlphaSimR. The clonal breeding framework was adapted from Bancic et al. (2025), *Plant breeding simulations with AlphaSimR*, Crop Science 65, e21312: https://doi.org/10.1002/csc2.21312. The original generic feature examples and wheat vignette have been removed from this almond-specific repository.

If you use this material, please cite the above publication and the original AlphaSimR publication:

    @article{Gaynor2021AlphaSimR,
      author = {Gaynor, R. Chris and Gorjanc, Gregor and Hickey, John M.},
      title = {{AlphaSimR}: An {R}-package for Breeding Program Simulations},
      volume = {11},
      year = {2021},
      doi = {10.1093/g3journal/jkaa017},
      journal = {G3: Genes|Genomes|Genetics}
    }



See [current experiment and external data storage](02_ClonalBreeding/STORAGE.md).

## Current manuscript design

The active catalogue contains 24 schemes with first parental use in years 9, 6, 4 or 2.
The manuscript comparison excludes the three exploratory GS4 OCS schemes (21 schemes).
See [the scenario catalogue](02_ClonalBreeding/SCHEMES.md) and
[external storage](02_ClonalBreeding/STORAGE.md). OCS uses optiSel contributions and mating
allocation with custom integer refinement, exact biological constraints and a maximum
combined genetic contribution of 10% per individual in new v8 runs.
Historical simulations belong to the previous design and are preserved externally.
The ongoing v7 experiment retains its frozen uncapped code. Updating this repository
does not start a simulation. Pending OCS in that experiment were cancelled by the user;
new capped OCS runs require a later explicit launch.

## Source layout

- `02_ClonalBreeding/main.R`: explicit entry point for a fresh complete catalogue run.
- `02_ClonalBreeding/00_Burn_in`: common trait/founder settings and standalone burn-in entry point.
- Numbered scenario folders: settings and explicit per-scenario launchers.
- `02_ClonalBreeding/redesigned_schemes.R`: annual simulation and decision-specific prediction models.
- `ocs_discrete.R`, `ocs_optisel.R`, `ocs_optisel_roles.R`, `ocs_cap10_support.R` and the compatibility adapter: active OCS backend.
- `plot_results.R`, `pedigree_representation.R` and `parallel_replicates.R`: plotting, pedigree export and parallel execution helpers.
- `tests/check_ocs_cap10.R`: small allocation checks; no breeding generations.

Generated results and retired local files remain outside Git. Source updates do not resume or launch any experiment.
