# Manuscript figures

`plot_manuscript.R` is a read-only plotting module. It never runs simulations,
fits breeding prediction models or computes new OCS allocations. `plot_results.R` exposes
it through the existing `plot_experiment()` entry point for future runs.

## Layout and scope

Each trajectory figure is a horizontal triptych: **selection accuracy**, **genetic gain**, **genetic variance**. Lines are solid.
Optional pointwise 95% Monte Carlo intervals require at least two replicates.

1. **Selection method:** PH9, PED9, GS9.
2. **Generation interval:** those three plus PH6, PED6, GS6, PED4, GS4 and SB2.
3. **Selection intensity:** GS6, PED6 and SB2, each with standard 500 and early
   1,500 -> 500 selection (suffix E).
4. **OCS:** GS6-E and SB2-E controls, each with Diversity (Ne=30), Balanced (Ne=10)
   and Gain (Ne=5).
5. **Final-year overview:** a vertical stack of two boxplots, genetic gain above genetic
   variance, with all 18 main-text schemes. Methods blocks determine order and fill:
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
are 18 x 6 inches; the endpoint grid is 18 x 9 inches. Outputs stay outside Git.

For an explicitly provisional example while the experiment is running:

```r
plot_manuscript("/path/to/run", preview = TRUE,
                out = "/path/to/external/preview")
```

Only full, unique 40-year trajectories enter the figures. Trajectory means use
shared complete replicate IDs within each panel. In preview mode, methods with
no results yet are excluded from pairing but remain in the legend.
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
legends and Methods-block labels remain. Sample sizes are exported in CSV files, not printed inside figures. Caption information,
provisional coverage and historical-design provenance are stored externally in
the exported README and CSV files rather than printed on the image.

## Legend columns

Trajectory legends read down each column. Selection-method and generation-interval
panels have separate Pheno, Pedigree and GS columns, with first-use years ordered
9, 6, 4, 2 where applicable. Empty alignment cells are invisible. Intensity uses
Pedigree, GS6 and speed-breeding columns, each showing the standard control above
early selection. OCS has GS6 and speed-breeding columns, each ordered control,
Diversity, Balanced, Gain. No overall figure titles or explanatory footnotes are added.

## OCS palette and endpoint metrics

The OCS trajectory panel uses blues for GS6 and turquoises for speed breeding,
while preserving control, Diversity, Balanced and Gain order in both columns.
Other trajectory panels retain the method-family palette. Actual panel colours
are exported in panel_color_key.csv. All trajectory panels read selection accuracy,
genetic gain, genetic variance from left to right. The endpoint grid contains only
genetic gain and genetic variance; prediction accuracy is omitted from that grid.


## Endpoint significance letters

The final-year gain and variance boxplots automatically receive compact letter
annotations when all 18 requested schemes have every planned paired replicate
(minimum three). Partial previews retain the boxes without letters;
the comparison family is never silently reduced to available schemes.

Each metric is analyzed separately with a REML mixed model: fixed scheme and
random replicate intercept (`nlme::lme`). Outcomes remain on their original
scale, with numerical rescaling only. The model assumes a common residual
variance and approximately normal conditional errors. Inspect the exported
normalized residuals and fitted models before publication; if these assumptions
are unsuitable, revise the analysis explicitly rather than treating the letters
as validated by default.

`emmeans` pairwise contrasts use Tukey correction across all main-text schemes
within each metric (alpha 0.05). `multcomp::cld`/`multcompView` provide letters;
sharing a letter means failure to detect a difference, not equivalence. Letters
are sorted with the highest estimated mean first; this is not a declaration
that higher variance is biologically preferable.

Outputs `endpoint_tests_{status,letters,pairs,means,diagnostics}.csv`,
`endpoint_tests_models.rds` and `endpoint_tests_method.txt` preserve the analysis.
Missing packages, incomplete coverage or fitting failures suppress the affected
metric's letters and are reported in the status file. Required packages: nlme,
emmeans, multcomp and multcompView. No simulations or OCS allocations are run.
