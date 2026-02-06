# pepdiff

Differential abundance analysis for phosphoproteomics data.

## Purpose

Identify peptides with significant abundance changes between
experimental conditions. Answers: “What’s differentially abundant?”

## Project Status

**v1.0.0 released** (Feb 2026). The package is feature-complete with: -
S3 classes: `pepdiff_data` and `pepdiff_results` - Three analysis
methods: GLM (Gamma + emmeans), ART, pairwise - Four pairwise tests:
wilcoxon, bootstrap_t, bayes_t, rankprod - GLM fit diagnostics with
[`plot_fit_diagnostics()`](https://teammaclean.github.io/pepdiff/reference/plot_fit_diagnostics.md) -
Plot methods for both classes - Six vignettes covering all workflows -
pkgdown documentation site: <https://teammaclean.github.io/pepdiff/> -
`devtools::check()` passes with 0 errors/warnings

**Companion to peppwR:** - peppwR: “How many samples do I need?” (power
analysis, planning) - pepdiff: “What’s differentially abundant?”
(analysis, results)

## Core Workflow

    CSV → read_pepdiff() → pepdiff_data → compare() → pepdiff_results → plots

## Documentation

- **Getting Started**: `vignettes/basic_workflow.Rmd`
- **GLM Analysis**: `vignettes/glm_analysis.Rmd`
- **ART Analysis**: `vignettes/art_analysis.Rmd`
- **Pairwise Tests**: `vignettes/pairwise_tests.Rmd`
- **Checking Model Fit**: `vignettes/checking_fit.Rmd`
- **Diagnostic Plots**: `vignettes/diagnostic_plots.Rmd`

Online: <https://teammaclean.github.io/pepdiff/>

## Design Constraints

### Cross-Sectional Factorial Designs Only

pepdiff assumes **independent biological replicates** in each cell of
the experimental design. NOT for longitudinal/repeated-measures studies
where the same biological unit is tracked over time.

- No random effects needed
- Simple GLM is appropriate
- Technical replicates must be combined explicitly before analysis

### Flexible Factors

Factor names are user-defined at import, not hardcoded. Common setups: -
treatment × timepoint - genotype × treatment - condition × timepoint ×
tissue

## S3 Classes

### `pepdiff_data`

Import/preprocessed data from
[`read_pepdiff()`](https://teammaclean.github.io/pepdiff/reference/read_pepdiff.md).

**Slots:** - `data` - tibble: peptide, gene_id, \[factors…\], bio_rep,
value - `factors` - character: names of factor columns - `design` -
tibble: unique factor combinations with n_reps, n_peptides -
`missingness` - tibble: peptide, na_rate, mnar_score, mean_abundance -
`peptides` - character vector - `call` - original function call

**Methods:** - [`print()`](https://rdrr.io/r/base/print.html) - Summary
of structure, missingness -
[`summary()`](https://rdrr.io/r/base/summary.html) - Detailed breakdown
by group - [`plot()`](https://rdrr.io/r/graphics/plot.default.html) -
Multi-panel: PCA + distributions + missingness

### `pepdiff_results`

Analysis results from
[`compare()`](https://teammaclean.github.io/pepdiff/reference/compare.md).
Results in **long format** (tidy).

**Slots:** - `results` - tibble (long): peptide, gene_id, comparison,
\[factor levels\], fold_change, log2_fc, test, p_value, fdr,
significant - `comparisons` - tibble: comparison definitions -
`method` - “glm”, “art”, or “pairwise” - `diagnostics` - tibble
(nested): peptide, converged, deviance, residuals, std_residuals,
fitted - `params` - list: alpha, fdr_method, formula, etc. - `data` -
the pepdiff_data object used - `call` - original function call

**Methods:** - [`print()`](https://rdrr.io/r/base/print.html) - Summary:
n comparisons, n significant -
[`summary()`](https://rdrr.io/r/base/summary.html) - Per-comparison
breakdown - [`plot()`](https://rdrr.io/r/graphics/plot.default.html) -
Multi-panel: volcano + p-value histogram + FC distribution

## Core Functions

### Import

``` r
read_pepdiff(file, id, gene, value, factors, replicate, tech_rep = NULL)
# Returns pepdiff_data

combine_tech_reps(data, fun = mean)
# Explicit step to average technical replicates
```

### Analysis

Three methods:

| Method            | Model                  | Use Case                   |
|-------------------|------------------------|----------------------------|
| `"glm"` (default) | Gamma GLM + emmeans    | Most proteomics data       |
| `"art"`           | Aligned Rank Transform | Non-parametric alternative |
| `"pairwise"`      | Direct two-group tests | Simple comparisons         |

**Simple interface:**

``` r
compare(data,
        compare = "treatment",       # factor to contrast
        ref = "ctrl",                # reference level
        within = "timepoint",        # stratify by (optional)
        method = "glm")
```

**Formula interface (power users):**

``` r
compare(data,
        contrast = pairwise ~ treatment | genotype + timepoint,
        method = "glm")
```

**Pairwise tests** (matching peppwR): - `"wilcoxon"` - Wilcoxon
rank-sum - `"bootstrap_t"` - Bootstrap t-test - `"bayes_t"` - Bayes
factor t-test - `"rankprod"` - Rank products

### Diagnostics

``` r
plot_fit_diagnostics(results)
# Returns 4-panel diagnostic plot for GLM model fit assessment
```

### Plots

**Class methods:** - `plot(pepdiff_data)` - PCA + distributions +
missingness - `plot(pepdiff_results)` - volcano + p-value hist + FC
distribution

**Individual functions:** -
[`plot_pca()`](https://teammaclean.github.io/pepdiff/reference/plot_pca.md),
`plot_distributions()`, `plot_missingness()` - `plot_volcano()`,
[`plot_heatmap()`](https://teammaclean.github.io/pepdiff/reference/plot_heatmap.md),
`plot_pvalue_hist()`, `plot_fc_distribution()` -
[`plot_fit_diagnostics()`](https://teammaclean.github.io/pepdiff/reference/plot_fit_diagnostics.md) -
GLM model fit assessment

## File Structure

    R/
      data.R            # read_pepdiff(), combine_tech_reps(), pepdiff_data class
      compare.R         # compare() generic and methods
      models.R          # GLM fitting, ART fitting, contrast extraction
      tests.R           # Pairwise statistical tests (wilcoxon, bootstrap_t, etc.)
      results.R         # pepdiff_results class, print/summary methods
      plots.R           # All plot functions and plot methods
      diagnostics.R     # plot_fit_diagnostics() and helpers
      utils.R           # Helpers, validation, internal utilities
      legacy.R          # Deprecated compare.data.frame method
      legacy-pepdiff.R  # Original v1 functions (preserved for compatibility)

    tests/testthat/
      helper-fixtures.R # Synthetic test data generators
      test-data.R       # Import and preprocessing tests
      test-compare.R    # compare() function tests
      test-models.R     # Model fitting tests
      test-tests.R      # Statistical test implementations
      test-results.R    # Results class tests
      test-plots.R      # Plot output tests
      test-diagnostics.R # Diagnostics function tests
      test-legacy.R     # Backwards compatibility tests

    vignettes/
      basic_workflow.Rmd   # Getting started guide
      glm_analysis.Rmd     # GLM method deep dive
      art_analysis.Rmd     # ART method guide
      pairwise_tests.Rmd   # Pairwise comparison methods
      checking_fit.Rmd     # GLM diagnostics guide
      diagnostic_plots.Rmd # Visualization options

## Error Handling

### Model Convergence Failures

If GLM/ART doesn’t converge for a peptide: - Peptide excluded from
results - Tracked in `diagnostics` slot - Warning in
[`print()`](https://rdrr.io/r/base/print.html): “X peptides excluded
(model did not converge)”

**Philosophy:** Fail is fail. User needs to know they may need a
different design.

### FDR Correction

Benjamini-Hochberg applied **within each comparison**, not globally
across all comparisons.

------------------------------------------------------------------------

## Development Guidelines

### Code Style

- Tidyverse style (pipes, dplyr verbs)
- S3 classes for core objects with print/plot methods
- Document with roxygen2
- Explicit namespace calls for non-base functions
  ([`dplyr::filter()`](https://dplyr.tidyverse.org/reference/filter.html))
- Keep test implementations in sync with peppwR (wilcoxon, bootstrap_t,
  bayes_t, rankprod)

### Package Checks

``` r
# Fast check (skip vignettes)
devtools::check(vignettes = FALSE)

# Run tests only
devtools::test()

# Full check
devtools::check()
```

## Dependencies

### Imports

- dplyr, tidyr, tibble, rlang, magrittr - data manipulation
- readr - CSV import
- ggplot2, cowplot - core plotting
- emmeans - GLM contrast extraction
- ARTool - ART method
- stringr, forcats - string/factor utilities

### Suggests

- ComplexHeatmap - heatmaps (Bioconductor)
- RankProd - rank products (Bioconductor)
- UpSetR - upset plots
- MKinfer - bootstrap tests
- testthat - testing
- knitr, rmarkdown - vignettes

Note: bayes_t uses a native JZS approximation (no BayesFactor
dependency)

## Development Workflow for Future Changes

See `semi-autonomous-feature-development.md` for detailed workflow.

### Discuss → TDD → Ralph Loop

1.  **Discuss** - Human describes feature/bug, Claude asks clarifying
    questions, agree on scope
2.  **TDD** - Write failing test first (the test IS the spec)
3.  **Ralph Loop** - Claude iterates autonomously until tests pass

### Key Principles

- **Tests are the contract** - No ambiguity about completion
- **No implementation until test fails** - Red → Green → Refactor
- **Clear context before implementation** - Commit test, start fresh
  session
- **Self-contained prompts** - Reference files, not discussion history

### For Bug Fixes / Features

    1. Discuss requirements
    2. Write failing test in tests/testthat/
    3. Commit the test
    4. /clear or new session
    5. /ralph-loop with verification command: devtools::test(filter = "test-name")
    6. Human smoke test
