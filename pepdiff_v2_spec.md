# pepdiff v2 Specification

Differential abundance analysis for PRM proteomics data.

## Purpose

Identify peptides with significant abundance changes between experimental conditions in cross-sectional factorial designs.

**Complements peppwR:**
- peppwR: "How many samples do I need?" (power analysis, planning)
- pepdiff: "What's differentially abundant?" (analysis, results)

## Design Philosophy

### Cross-Sectional Factorial Designs Only

pepdiff assumes **independent biological replicates** in each cell of the experimental design. It is **not** intended for longitudinal/repeated-measures studies where the same biological unit is tracked over time.

Factor names are user-defined. Common setups:
- treatment × timepoint
- genotype × treatment
- condition × timepoint × tissue

Technical replicates (same sample measured multiple times) should be combined using `combine_tech_reps()` before analysis.

### Analysis Approach

Three analysis methods:

| Method | Model | Contrasts | Use Case |
|--------|-------|-----------|----------|
| `"glm"` (default) | Gamma GLM with log link | emmeans | Most proteomics data (positive, right-skewed) |
| `"art"` | Aligned Rank Transform | art.con() | Non-parametric alternative, distribution-free |
| `"pairwise"` | None (direct tests) | N/A | Legacy mode, simple two-group comparisons |

GLM coefficients on the log scale are directly interpretable as log fold changes.

### Tidy Data Principles

Results stored in long format (one row per peptide × comparison × test), not wide format with columns like `wilcoxon_p_val`, `bootstrap_t_p_val`. This enables easy filtering, summarizing, and plotting with tidyverse tools.

---

## S3 Classes

### `pepdiff_data`

Created by `read_pepdiff()`. Contains imported and preprocessed data.

**Slots:**
```r
├── data            # tibble: peptide, gene_id, [factors...], bio_rep, value
├── factors         # character: names of factor columns
├── design          # tibble: unique factor combinations with n_reps, n_peptides
├── missingness     # tibble: peptide, na_rate, mnar_score, mean_abundance
├── peptides        # character: unique peptide IDs
└── call            # original function call
```

**Methods:**
- `print()` - Summary: n peptides, factors, replicates per group, missingness rate
- `summary()` - Detailed breakdown by group, missingness stats, MNAR warnings
- `plot()` - Multi-panel diagnostic: PCA + distributions + missingness heatmap

### `pepdiff_results`

Created by `compare()`. Contains analysis results in long format.

**Slots:**
```r
├── results         # tibble (long format):
│                   #   peptide, gene_id, comparison, [factor levels],
│                   #   fold_change, log2_fc, test, p_value, fdr, significant
├── comparisons     # tibble: comparison definitions that were run
├── method          # character: "glm", "art", or "pairwise"
├── diagnostics     # tibble (nested): peptide, converged, deviance, model, residuals
├── params          # list: alpha, fdr_method, formula used, etc.
├── data            # the pepdiff_data object used
└── call            # original function call
```

**Methods:**
- `print()` - Summary: n comparisons, n significant peptides (at default threshold)
- `summary()` - Per-comparison, per-test breakdown of significant counts
- `plot()` - Multi-panel summary: volcano + p-value histogram + FC distribution

---

## Core Functions

### Import & Preprocessing

#### `read_pepdiff()`

Import data from CSV file.

```r
read_pepdiff(
  file,
  id = "peptide",
  gene = "gene_id",
  value = "abundance",
  factors,                    # required: character vector of factor column names
  replicate = "bio_rep",
  tech_rep = NULL             # if present, warns user to combine first
)
```

**Returns:** `pepdiff_data` object

**Behavior:**
- Validates required columns exist
- Computes missingness statistics per peptide (NA rate, MNAR score)
- Builds design summary (factor combinations × replicate counts)
- If `tech_rep` specified and multiple tech reps detected, warns user to run `combine_tech_reps()`

#### `combine_tech_reps()`

Average technical replicates.
```r
combine_tech_reps(data, fun = mean)
```

**Arguments:**
- `data` - `pepdiff_data` object with tech_rep column
- `fun` - aggregation function (default: `mean`)

**Returns:** `pepdiff_data` object with tech reps collapsed

---

### Analysis

#### `compare()`
Main analysis function. Generic with methods for different input types.

**Simple interface** (common cases):
```r
compare(
  data,                       # pepdiff_data object
  compare,                    # character: factor to contrast
  ref,                        # character: reference level
  within = NULL,              # character vector: factors to stratify by
  method = "glm",             # "glm", "art", or "pairwise"
  test = "wilcoxon",          # for pairwise method only
  alpha = 0.05,
  fdr_method = "BH"
)
```

**Formula interface** (power users):
```r
compare(
  data,
  contrast,                   # formula: e.g., pairwise ~ treatment | timepoint
  method = "glm",
  alpha = 0.05,
  fdr_method = "BH"
)
```

**Returns:** `pepdiff_results` object

**Examples:**
```r
# Simple: treatment vs control, within each timepoint
compare(dat, compare = "treatment", ref = "ctrl", within = "timepoint")

# Two factors to stratify by
compare(dat, compare = "treatment", ref = "ctrl", within = c("genotype", "timepoint"))

# Formula interface
compare(dat, contrast = pairwise ~ treatment | genotype + timepoint)

# Non-parametric (ART)
compare(dat, compare = "treatment", ref = "ctrl", method = "art")

# Legacy pairwise tests
compare(dat, compare = "treatment", ref = "ctrl", method = "pairwise", test = "wilcoxon")
compare(dat, compare = "treatment", ref = "ctrl", method = "pairwise", test = c("wilcoxon", "bootstrap_t"))
```

### Statistical Tests (pairwise method)

Matching peppwR for consistency:

| Test | ID | Description |
|------|----|-------------|
| Wilcoxon rank-sum | `"wilcoxon"` | Non-parametric, robust |
| Bootstrap-t | `"bootstrap_t"` | Resampling-based, handles non-normality |
| Bayes factor t-test | `"bayes_t"` | Bayesian evidence for effect |
| Rank Products | `"rankprod"` | Designed for omics, handles small samples |

---

## Plots

### Class Methods (Multi-Panel Summaries)

#### `plot.pepdiff_data()`
Exploratory data visualization.

**Panels:**
1. PCA of samples (colored by factors)
2. Abundance distributions per group
3. Missingness heatmap

#### `plot.pepdiff_results()`
Results summary visualization.

**Panels:**
1. Volcano plot (log2 FC vs -log10 p-value)
2. P-value histogram
3. Fold change distribution

### Individual Plot Functions

For fine-grained control and single-purpose plots.

**Data exploration:**
```r
plot_pca(data, color_by = NULL)
plot_distributions(data, facet_by = NULL)
plot_missingness(data)
```

**Results visualization:**
```r
plot_volcano(results, comparison = NULL, alpha = 0.05, fc_threshold = 1)
plot_heatmap(results, cluster_rows = TRUE, cluster_cols = TRUE)
plot_upset(results, min_size = 1)
plot_pvalue_hist(results, comparison = NULL)
plot_fc_distribution(results, comparison = NULL)
```

---

## Error Handling

### Model Convergence Failures

For GLM and ART methods, some peptides may fail to converge (insufficient data, separation, etc.).

**Behavior:**
- Failed peptides excluded from results
- Tracked in `diagnostics` slot: `results$diagnostics |> filter(!converged)`
- Clear warning in `print()`: "⚠ 23 peptides excluded (model did not converge)"

**Philosophy:** Fail is fail. If a peptide can't be modeled with the given design, the user needs to know, not receive imputed/fallback results.

### Multiple Testing Correction

FDR correction (Benjamini-Hochberg by default) applied **within each comparison**, not across all comparisons globally.

---

## File Structure

```
R/
  data.R            # read_pepdiff(), combine_tech_reps(), pepdiff_data class
  compare.R         # compare() generic and methods
  models.R          # GLM fitting, ART fitting, contrast extraction
  tests.R           # Pairwise statistical tests (wilcoxon, bootstrap_t, etc.)
  results.R         # pepdiff_results class, print/summary methods
  plots.R           # All plot functions and plot methods
  utils.R           # Helpers, validation, internal utilities

tests/testthat/
  test-data.R       # Import and preprocessing tests
  test-compare.R    # compare() function tests
  test-models.R     # Model fitting tests
  test-tests.R      # Statistical test implementations
  test-results.R    # Results class tests
  test-plots.R      # Plot output tests

vignettes/
  getting-started.Rmd      # Basic workflow introduction
  analysis-workflow.Rmd    # Complete worked example
  advanced-designs.Rmd     # Multi-factor designs, formula interface
```

---

## Dependencies

### Imports
- dplyr, tidyr, tibble, purrr, rlang - data manipulation
- readr - CSV import
- ggplot2, cowplot - core plotting
- emmeans - GLM contrast extraction
- stats - base statistical functions

### Suggests
- ARTool - Aligned Rank Transform method
- BayesFactor - Bayes factor t-test (pairwise method)
- ComplexHeatmap - heatmap plots
- UpSetR - UpSet plots
- RColorBrewer, viridis - color palettes
- testthat - testing
- knitr, rmarkdown - vignettes

---

## Development Workflow

### TDD + Ralph Loop

1. **Discuss** - Clarify function purpose, interface, edge cases
2. **Specify** - Write failing tests that define expected behavior
3. **Implement** - Ralph Loop until tests pass
4. **Document** - Roxygen comments, update vignettes if needed

### Test Coverage

Every exported function must have tests covering:
- Happy path (typical usage)
- Edge cases (empty input, single peptide, single replicate)
- Error conditions (invalid input, missing columns)
- Return value structure validation

---

## Version History

- **v0.1.x** - Prototype (current pepdiff-master)
- **v2.0.0** - Complete rewrite per this specification

---

## Relationship to peppwR

These packages share:
- Statistical test implementations (wilcoxon, bootstrap_t, bayes_t, rankprod)
- Missingness detection approach (NA rate, MNAR score)
- Code style (tidyverse, explicit namespacing, S3 classes)
- Documentation approach (pedagogic vignettes)

They are **separate repositories** but should be kept in sync manually. Changes to shared test implementations should be mirrored.