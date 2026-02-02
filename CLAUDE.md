# pepdiff

Differential abundance analysis for PRM proteomics data.

## Purpose

Identify peptides with significant abundance changes between experimental conditions. Answers: "What's differentially abundant?"

**Companion to peppwR:**
- peppwR: "How many samples do I need?" (power analysis, planning)
- pepdiff: "What's differentially abundant?" (analysis, results)

## Core Workflow

```
CSV → read_pepdiff() → pepdiff_data → compare() → pepdiff_results → plots
```

## Design Constraints

### Cross-Sectional Factorial Designs Only

pepdiff assumes **independent biological replicates** in each cell of the experimental design. NOT for longitudinal/repeated-measures studies where the same biological unit is tracked over time.

- No random effects needed
- Simple GLM is appropriate
- Technical replicates must be combined explicitly before analysis

### Flexible Factors

Factor names are user-defined at import, not hardcoded. Common setups:
- treatment × timepoint
- genotype × treatment
- condition × timepoint × tissue

## S3 Classes

### `pepdiff_data`

Import/preprocessed data from `read_pepdiff()`.

**Slots:**
- `data` - tibble: peptide, gene_id, [factors...], bio_rep, value
- `factors` - character: names of factor columns
- `design` - tibble: unique factor combinations with n_reps, n_peptides
- `missingness` - tibble: peptide, na_rate, mnar_score, mean_abundance
- `peptides` - character vector
- `call` - original function call

**Methods:**
- `print()` - Summary of structure, missingness
- `summary()` - Detailed breakdown by group
- `plot()` - Multi-panel: PCA + distributions + missingness

### `pepdiff_results`

Analysis results from `compare()`. Results in **long format** (tidy).

**Slots:**
- `results` - tibble (long): peptide, gene_id, comparison, [factor levels], fold_change, log2_fc, test, p_value, fdr, significant
- `comparisons` - tibble: comparison definitions
- `method` - "glm", "art", or "pairwise"
- `diagnostics` - tibble (nested): peptide, converged, deviance, model, residuals
- `params` - list: alpha, fdr_method, formula, etc.
- `data` - the pepdiff_data object used
- `call` - original function call

**Methods:**
- `print()` - Summary: n comparisons, n significant
- `summary()` - Per-comparison breakdown
- `plot()` - Multi-panel: volcano + p-value histogram + FC distribution

## Core Functions

### Import

```r
read_pepdiff(file, id, gene, value, factors, replicate, tech_rep = NULL)
# Returns pepdiff_data

combine_tech_reps(data, fun = mean)
# Explicit step to average technical replicates
```

### Analysis

Three methods:

| Method | Model | Use Case |
|--------|-------|----------|
| `"glm"` (default) | Gamma GLM + emmeans | Most proteomics data |
| `"art"` | Aligned Rank Transform | Non-parametric alternative |
| `"pairwise"` | Direct two-group tests | Legacy mode |

**Simple interface:**
```r
compare(data,
        compare = "treatment",       # factor to contrast
        ref = "ctrl",                # reference level
        within = "timepoint",        # stratify by (optional)
        method = "glm")
```

**Formula interface (power users):**
```r
compare(data,
        contrast = pairwise ~ treatment | genotype + timepoint,
        method = "glm")
```

**Pairwise tests** (matching peppwR):
- `"wilcoxon"` - Wilcoxon rank-sum
- `"bootstrap_t"` - Bootstrap t-test
- `"bayes_t"` - Bayes factor t-test
- `"rankprod"` - Rank products

### Plots

**Class methods:**
- `plot(pepdiff_data)` - PCA + distributions + missingness
- `plot(pepdiff_results)` - volcano + p-value hist + FC distribution

**Individual functions:**
- `plot_pca()`, `plot_distributions()`, `plot_missingness()`
- `plot_volcano()`, `plot_heatmap()`, `plot_upset()`, `plot_pvalue_hist()`, `plot_fc_distribution()`

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
```

## Error Handling

### Model Convergence Failures

If GLM/ART doesn't converge for a peptide:
- Peptide excluded from results
- Tracked in `diagnostics` slot
- Warning in `print()`: "⚠ X peptides excluded (model did not converge)"

**Philosophy:** Fail is fail. User needs to know they may need a different design.

### FDR Correction

Benjamini-Hochberg applied **within each comparison**, not globally across all comparisons.

---

## Development Guidelines

### Code Style

- Tidyverse style (pipes, dplyr verbs)
- S3 classes for core objects with print/plot methods
- Document with roxygen2
- Explicit namespace calls for non-base functions (`dplyr::filter()`)
- Keep test implementations in sync with peppwR (wilcoxon, bootstrap_t, bayes_t, rankprod)

### Development Workflow

We use the **Discuss → TDD → Ralph Loop** workflow. See `semi-autonomous-feature-development.md` for details.

**Summary:**
1. **Discuss**: Reach shared understanding of the task/feature
2. **TDD**: Write failing test that captures the requirement, commit it
3. **Clear Context**: `/clear` or new session to maximize implementation context
4. **Ralph Loop**: Autonomous implementation with self-contained prompt
5. **Smoke Test**: Human final verification

### Ralph Loop Prompt Template

```
/ralph-loop "Implement [feature] for pepdiff.

## Failing Test
tests/testthat/test-[feature].R

## Relevant Files
- R/[file].R
- [other files]

## Verification
Rscript -e 'devtools::test(filter = "[feature]")'
Must show: OK

## Success Criteria
All tests pass, devtools::check() has 0 errors." --completion-promise "FEATURE-COMPLETE" --max-iterations 20
```

### Key Principles

- **Tests are the contract**: No ambiguity about completion
- **Small tasks**: Better speed and accuracy than large tasks
- **Context in files**: Specs live in test files, not conversation history
- **Clear before execute**: Maximize context for implementation work

## Package Checks

### Routine Development (fast)

Skip vignettes during routine checks:

```r
devtools::check(vignettes = FALSE)
```

Or from command line:

```bash
Rscript -e "devtools::check(vignettes = FALSE)"
```

### Run Tests Only

```r
devtools::test()
# Or specific test file:
devtools::test(filter = "compare")
```

### Full Check (before release)

```r
devtools::build_vignettes()
devtools::check()
```

## Dependencies

### Imports
- dplyr, tidyr, tibble, purrr, rlang - data manipulation
- readr - CSV import
- ggplot2, cowplot - core plotting
- emmeans - GLM contrast extraction

### Suggests
- ARTool - ART method
- BayesFactor - Bayes factor t-test
- ComplexHeatmap, UpSetR - specialized plots
- testthat - testing
- knitr, rmarkdown - vignettes

---

## Reference

Full specification: `pepdiff_v2_spec.md`