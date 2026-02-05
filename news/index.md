# Changelog

## pepdiff 1.0.0

Initial release of pepdiff for differential abundance analysis of PRM
proteomics data.

### Features

#### Analysis Methods

- [`compare()`](https://teammaclean.github.io/pepdiff/reference/compare.md) -
  Main analysis function supporting three methods:
  - **GLM**: Gamma GLM with emmeans for factorial designs (recommended)
  - **ART**: Aligned Rank Transform for non-parametric analysis
  - **Pairwise**: Direct two-group comparisons (Wilcoxon, bootstrap-t,
    Bayes factor, rank products)
- Stratified comparisons with `within` parameter for analysing effects
  within factor levels

#### Data Management

- [`read_pepdiff()`](https://teammaclean.github.io/pepdiff/reference/read_pepdiff.md) -
  Import CSV data with flexible factor specification
- [`combine_tech_reps()`](https://teammaclean.github.io/pepdiff/reference/combine_tech_reps.md) -
  Combine technical replicates before analysis
- S3 classes `pepdiff_data` and `pepdiff_results` with print, summary,
  and plot methods

#### Diagnostics

- [`plot_fit_diagnostics()`](https://teammaclean.github.io/pepdiff/reference/plot_fit_diagnostics.md) -
  Four-panel diagnostic plot for assessing GLM model fit
- Stored residuals and fitted values for post-hoc diagnostics
- Convergence tracking for all fitted models

#### Visualization

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  both data and results objects
- `plot_volcano()` - Volcano plots with customizable thresholds
- [`plot_heatmap()`](https://teammaclean.github.io/pepdiff/reference/plot_heatmap.md) -
  Heatmaps of significant peptides (requires ComplexHeatmap)
- [`plot_pca()`](https://teammaclean.github.io/pepdiff/reference/plot_pca.md) -
  PCA visualization of samples
- `plot_pvalue_hist()` - P-value distribution histograms
- `plot_fc_distribution()` - Fold change distributions
- `plot_missingness()` - Missing data patterns
- `plot_distributions()` - Abundance distributions by group

#### Companion Package

pepdiff is designed to work alongside
[peppwR](https://github.com/TeamMacLean/peppwR): - **peppwR**: “How many
samples do I need?” (power analysis) - **pepdiff**: “What’s
differentially abundant?” (analysis)
