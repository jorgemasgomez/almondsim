# Manuscript figures

`plot_manuscript.R` is a read-only plotting module. It never runs simulations,
fits prediction models or computes new OCS allocations. `plot_results.R` exposes
it through the existing `plot_experiment()` entry point for future runs.

## Layout and scope

Each trajectory figure is a horizontal triptych: **genetic gain**, **genetic
variance (measured in seedlings)**, **selection accuracy**. Lines are solid.
Optional pointwise 95% Monte Carlo intervals require at least two replicates.

1. **Selection method:** PH9, PED9, GS9.
2. **Generation interval:** those three plus PH6, PED6, GS6, PED4, GS4 and SB2.
3. **Selection intensity:** GS6, PED6 and SB2, each with standard 500 and early
   1,500 -> 500 selection (suffix E).
4. **OCS:** GS6-E and SB2-E controls, each with Diversity (.10), Balanced (.20)
   and Gain (.40). Exploratory GS4 OCS is excluded.
5. **Supplement S01:** GS6 renewal 25%, 50%, 75% and 100%. These variants are
   kept out of the main trajectory panels and the final-year overview.
6. **Final-year overview:** a vertical stack of three boxplots, in the same metric
   order, with all 18 main-text schemes. Methods blocks determine order and fill:
   selection accuracy, generation interval, selection intensity, OCS.

Trajectory colours are fixed by method family: Pheno orange, pedigree purple,
GS blue. Variants use different shades. Boxplot colours instead identify Methods
blocks. Abbreviations and per-scheme colours are exported in `scheme_key.csv`.
A dot in the endpoint grid represents one replicate. Single-replicate schemes
show that observation; unavailable schemes remain empty. PH accuracy is NA.

## Use

```r
source("02_ClonalBreeding/plot_manuscript.R")
plot_manuscript("/path/to/run", with_intervals = TRUE)
plot_manuscript("/path/to/run", with_intervals = FALSE)
```

The default output directory is `run/manuscript_figures/with_CI` or `lines`.
PNG files are 400 dpi; PDF files retain vector plot elements. Trajectory figures
are 18 x 6 inches; the endpoint grid is 18 x 12.5 inches. Outputs stay outside Git.

For an explicitly provisional example while the experiment is running:

```r
plot_manuscript("/path/to/run", preview = TRUE,
                out = "/path/to/external/preview")
```

Only full, unique 40-year trajectories enter the figures. Trajectory means use
shared complete replicate IDs within each panel. In preview mode, methods with
no results yet are excluded from pairing but remain in the legend with n=0.
In final mode, all requested methods must have a shared complete replicate.
`paired=FALSE` optionally displays per-scheme complete-replicate means, with
unequal coverage explicitly identified. Endpoint boxplots show each scheme's
available complete replicates, not paired differences.

Coverage, plotted summaries, individual endpoint values and figure data are
exported alongside the images. Gain is relative to each replicate's burn-in.
Seedling gain and variance are post-early-selection; prediction accuracy is
pre-early-selection. No zero accuracy is invented for phenotypic selection.

The run's recorded design identifier is respected. Historical OCS results carry
an explicit historical-design note and are never relabelled as native-offspring
OCS. Data from different runs or design versions are not combined.
Existing frozen simulation code is not modified by changing this plot module.

## Publication styling

Figures have no overall title, subtitle or explanatory footnote. Metric labels
are Genetic gain, Genetic variance and Selection accuracy. Axis labels, scheme
legends, Methods-block labels and sample sizes remain. Caption information,
provisional coverage and historical-design provenance are stored externally in
the exported README and CSV files rather than printed on the image.

## Legend columns

Trajectory legends read down each column. Selection-method and generation-interval
panels have separate Pheno, Pedigree and GS columns, with first-use years ordered
9, 6, 4, 2 where applicable. Empty alignment cells are invisible. Intensity uses
Pedigree, GS6 and speed-breeding columns, each showing the standard control above
early selection. OCS has GS6 and speed-breeding columns, each ordered control,
Diversity, Balanced, Gain. The supplementary renewal legend reads 25%, 50%, 75%,
100%. No overall figure titles or explanatory footnotes are added.
