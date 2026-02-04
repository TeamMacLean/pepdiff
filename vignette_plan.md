# pepdiff Vignette Plan

## Audience

Intelligent but statistically uninformed. Familiar with hypothesis testing basics (p-values, significance) and understands simple linear models conceptually. Wants practical guidance, not theory lectures.

## Data Strategy

- **Vignettes**: Synthetic data generated within the vignette (reproducible, no external dependencies)
- **Tutorials** (future): Pre-rendered with larger real datasets, hosted externally

## Philosophy

- **Model fit over method selection**: Don't prescribe "use X when Y". Instead: fit the model, check diagnostics, if it fits then it's appropriate.
- **No imputation**: Accept that missing data means some peptides can't be analyzed. Refer to peppwR for understanding missingness patterns and power implications.
- **Fail gracefully**: When there aren't enough replicates, the analysis fails - that's information, not a bug.

---

## Vignette 1: `basic_workflow.Rmd`

**Title**: Getting Started with pepdiff

**Goal**: Complete workflow from import to results in ~10 minutes reading time

### Outline

1. **What pepdiff does**
   - One paragraph: differential abundance analysis for PRM proteomics
   - Companion to peppwR (planning) → pepdiff (analysis)

2. **Installation**
   - Standard install instructions
   - Note on Bioconductor dependencies (ComplexHeatmap, etc.)

3. **Importing data**
   - `read_pepdiff()` with synthetic data
   - Explain the column mapping: id, gene, value, factors, replicate
   - Show the `pepdiff_data` object structure
   - `print()` and `summary()` output

4. **Quick look at your data**
   - `plot(data)` - the multi-panel diagnostic
   - Brief interpretation: "PCA should show groups separating", "distributions should be similar across reps"

5. **Running a comparison**
   - Simple `compare(data, compare = "treatment", ref = "ctrl")`
   - Default is GLM - that's fine for most cases
   - Show the `pepdiff_results` object

6. **Understanding results**
   - Key columns: peptide, fold_change, p_value, fdr, significant
   - `significant()` accessor
   - `plot(results)` - volcano, p-value histogram

7. **Next steps**
   - Link to other vignettes for deeper dives

### Synthetic Data

```r
# 2 treatments, 3 reps, 50 peptides, ~30% differentially abundant
set.seed(123)
# Generate and write to temp CSV within vignette
```

---

## Vignette 2: `pairwise_tests.Rmd`

**Title**: Two-Group Comparisons with Pairwise Tests

**Goal**: When and how to use the pairwise method; understand the four test options

### Outline

1. **When to use pairwise tests**
   - Simple two-group comparison (ctrl vs treatment)
   - No complex factorial structure
   - Often sufficient for straightforward experiments

2. **The four tests**

   Each gets a subsection with:
   - One-sentence description
   - When it works well
   - Example with same data

   **Wilcoxon rank-sum** (`test = "wilcoxon"`)
   - Non-parametric, compares ranks not values
   - Robust to outliers and non-normality
   - Default choice for small samples

   **Bootstrap t-test** (`test = "bootstrap_t"`)
   - Resampling-based, doesn't assume normal distribution
   - Good when you have 4+ replicates per group
   - More powerful than Wilcoxon when data is roughly symmetric

   **Bayes factor t-test** (`test = "bayes_t"`)
   - Returns evidence for/against difference (not just "reject/don't reject")
   - BF > 3: moderate evidence, BF > 10: strong evidence
   - Useful when you want to quantify evidence, not just p-values

   **Rank products** (`test = "rankprod"`)
   - Specifically designed for high-throughput data
   - Handles asymmetric changes well
   - Returns separate p-values for up/down regulation

3. **Comparing results across tests**
   - Run all four on same data
   - Table of results
   - They often agree; when they don't, it's informative

4. **FDR correction**
   - What it is (controlling false discovery rate)
   - Applied within each comparison
   - Use `fdr` column, not raw `p_value`, for significance calls

5. **When pairwise isn't enough**
   - Multiple factors → consider GLM
   - Want to model interactions → need GLM

### Synthetic Data

```r
# Clear ctrl vs trt difference, 4 reps each, 30 peptides
# Some peptides strongly different, some moderately, some not at all
```

---

## Vignette 3: `glm_analysis.Rmd`

**Title**: Factorial Designs with GLM

**Goal**: Using Gamma GLM for experiments with multiple factors

### Outline

1. **When you need GLM**
   - Multiple experimental factors (treatment × timepoint × genotype)
   - Want to account for interactions
   - More replicates available (GLM needs more data)

2. **Why Gamma GLM?**
   - Proteomics abundances are positive, right-skewed
   - Gamma distribution models this naturally
   - Log link gives multiplicative (fold change) interpretation
   - Brief note: zeros will cause failure (that's correct behavior)

3. **Setting up a factorial analysis**
   - Import with multiple factors
   - `compare(data, compare = "treatment", ref = "ctrl", method = "glm")`
   - The `within` argument for stratified comparisons

4. **Understanding the output**
   - Contrasts from emmeans
   - Fold changes are back-transformed from log scale
   - Diagnostics slot: which peptides converged?

5. **Handling convergence failures**
   - Some peptides won't converge - that's OK
   - Common reasons: too few observations, zero variance, extreme values
   - Check `diagnostics` slot
   - Don't force it - if it fails, you don't have enough data for that peptide
   - Refer to peppwR for power planning

6. **Checking model fit**
   - Are residuals reasonable?
   - Does the volcano plot look sensible?
   - If many peptides fail, reconsider experimental design

7. **Interaction effects** (expanded)
   - What the warning means: "Results may be misleading due to involvement in interactions"
     - Main effect estimates are averaged across levels of other factors
     - If interaction exists, the average may not represent any real condition
   - When to worry:
     - Interaction term is significant
     - Effect direction differs across levels of another factor
     - Example: treatment helps at 0h but hurts at 24h → main effect averages to ~nothing
   - What to do:
     - Look at simple effects (treatment effect *within* each timepoint)
     - Use `within` argument: `compare(data, compare = "treatment", ref = "ctrl", within = "timepoint")`
     - Report the conditional effects, not the marginal main effect
   - When it's fine to ignore:
     - Interaction is small/non-significant
     - Effect is consistent in direction across all levels

### Synthetic Data

```r
# 2 treatments × 2 timepoints, 3 reps per cell, 40 peptides
# Treatment effect, timepoint effect, some interaction
```

---

## Vignette 4: `art_analysis.Rmd`

**Title**: Non-parametric Factorial Analysis with ART

**Goal**: Using Aligned Rank Transform when GLM assumptions are suspect

### Outline

1. **What is ART?**
   - Aligned Rank Transform: non-parametric method for factorial designs
   - Handles non-normal data without transformation
   - Maintains ability to test interactions

2. **When to consider ART**
   - GLM diagnostics look poor
   - Data is ordinal or highly non-normal
   - You're uncomfortable with distributional assumptions
   - Note: requires ARTool package

3. **Running ART analysis**
   - `compare(data, compare = "treatment", ref = "ctrl", method = "art")`
   - Same interface as GLM

4. **ART vs GLM comparison**
   - Run both on same data
   - Compare results
   - Usually similar when GLM is appropriate
   - ART may be more robust when GLM struggles

5. **Limitations**
   - Slower than GLM
   - Interpretation less intuitive (no direct fold changes)
   - Requires balanced-ish designs

6. **Practical recommendation**
   - Start with GLM (default)
   - If diagnostics suggest problems, try ART
   - If both fail, you may need more data

### Synthetic Data

```r
# Same factorial structure as GLM vignette
# But with some non-normal noise patterns
```

---

## Vignette 5: `diagnostic_plots.Rmd`

**Title**: Visual Quality Control and Interpretation

**Goal**: Use plots to assess data quality and interpret results

### Outline

1. **Why visualize?**
   - Numbers lie; plots reveal
   - QC before analysis, interpretation after

2. **Data-level diagnostics** (`plot(pepdiff_data)`)

   **PCA plot**
   - What to look for: replicates cluster, groups separate
   - Warning signs: outlier samples, no group separation, batch effects
   - What to do: investigate outliers, consider excluding bad samples

   **Distribution plots**
   - What to look for: similar shapes across samples
   - Warning signs: shifted distributions, different spreads
   - What to do: check for normalization issues, sample prep problems

   **Missingness plot**
   - What to look for: random scatter of missing values
   - Warning signs: systematic missingness (MNAR), one sample missing everything
   - What to do: refer to peppwR for missingness patterns; accept that high-missingness peptides may not be analyzable

3. **Results-level diagnostics** (`plot(pepdiff_results)`)

   **Volcano plot**
   - What to look for: symmetric spread, significant hits at edges
   - Warning signs: all hits in one direction, strange patterns
   - Customization: `fc_threshold`, `alpha`

   **P-value histogram**
   - What to look for: uniform with spike near 0 (if true positives exist)
   - Warning signs: U-shape (p-value inflation), spike at 1 (conservative)
   - What it means for your FDR estimates

   **Fold change distribution**
   - What to look for: centered near 0 (log2 scale)
   - Warning signs: systematic shift, bimodal distribution

4. **Individual plot functions**
   - `plot_pca_simple()`, `plot_volcano_new()`, etc.
   - When you want just one plot, not the grid
   - Customization options

5. **Exporting for publication**
   - These are ggplot2 objects
   - `ggsave()` workflow
   - Customizing themes

### Synthetic Data

```r
# Create "good" dataset and "problematic" dataset
# Show what each looks like in the plots
```

---

## Implementation Notes

### Order of writing
1. `basic_workflow` - foundation, others reference it
2. `pairwise_tests` - simplest method, natural second step
3. `diagnostic_plots` - useful for all methods
4. `glm_analysis` - more advanced
5. `art_analysis` - most specialized

### Common elements
- Each vignette starts with `library(pepdiff)`
- Use consistent synthetic data generation (can share helper functions)
- End with "See also" links to related vignettes
- Keep code chunks short and commented

### Tone
- Direct and practical: "Here's what to do"
- But also pedagogical: "Here's why it works"
- Not formal theory, but intuition that helps users understand and adapt
- Acknowledge uncertainty: "if this doesn't look right, investigate"
- No false precision: methods are tools, not truth machines

### Pedagogical notes to weave in
- **Gamma GLM**: Why right-skewed? Multiplicative processes, can't go negative, variance scales with mean. Log link → fold changes.
- **FDR**: If you call 100 things significant at 5% FDR, expect ~5 false positives. That's the deal you're making.
- **P-value histogram**: Uniform = null is true (no signal). Spike near 0 = true positives pulling p-values down. U-shape = something's wrong.
- **Wilcoxon**: Compares ranks, not values. Robust because outliers can only affect rank by 1.
- **Bootstrap**: Resamples your actual data to build a null distribution. Doesn't assume normality because it uses *your* distribution.
- **Bayes factor**: Quantifies evidence (how much more likely is H1 than H0?) rather than just reject/don't reject.
- **ART**: Aligns data to remove main effects, then ranks. Lets you test interactions non-parametrically.
- **Interactions**: Main effect = average across levels. If effect flips direction, the average is meaningless.

### Testing vignettes
- Build with `devtools::build_vignettes()`
- Ensure all code runs without external data
- Keep execution time reasonable (<30s per vignette)
