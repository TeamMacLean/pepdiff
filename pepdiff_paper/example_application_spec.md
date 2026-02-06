# pepdiff Example Application Specification

This document specifies the simulated dataset and analyses for the pepdiff paper's Example Application section and Figure 1.

## Overview

We demonstrate pepdiff using a simulated phosphoproteomics experiment with:
- Factorial design: treatment (ctrl, drug) × timepoint (0h, 6h, 24h)
- Realistic data characteristics (Gamma-distributed, heteroscedastic, MNAR)
- Known ground truth for validation

---

## 1. Simulated Dataset

### 1.1 Experimental Design

| Factor | Levels | Description |
|--------|--------|-------------|
| treatment | ctrl, drug | Two-group comparison |
| timepoint | 0h, 6h, 24h | Temporal dynamics |
| bio_rep | 1-4 | Biological replicates per condition |

**Total samples**: 2 × 3 × 4 = 24 samples

### 1.2 Peptide Parameters

Generate 500 peptides with heterogeneous parameters:

```r
set.seed(2026)

n_peptides <- 500
n_true_positives <- 50  # Peptides with real treatment effect

peptide_params <- tibble(
  peptide = paste0("pep_", sprintf("%04d", 1:n_peptides)),
  gene_id = paste0("gene_", sprintf("%03d", rep(1:100, each = 5))),


  # Gamma distribution parameters (heterogeneous)
  shape = runif(n_peptides, 1.5, 6),
  rate = runif(n_peptides, 0.03, 0.15),

  # True effect (only for first 50 peptides, only at 24h)
  true_effect = c(rep(TRUE, n_true_positives), rep(FALSE, n_peptides - n_true_positives)),
  effect_size = ifelse(true_effect, runif(n_peptides, 1.8, 2.5), 1.0),

  # Missingness probability (correlated with abundance = shape/rate)
  mean_abundance = shape / rate,
  miss_prob = pmax(0, 0.25 - 0.005 * mean_abundance)  # MNAR: low abundance = more missing
)
```

### 1.3 Data Generation

For each peptide, sample, and condition:

```r
generate_sample <- function(shape, rate, effect, miss_prob) {
  # Draw from Gamma
  value <- rgamma(1, shape = shape, rate = rate)

  # Apply effect (multiplicative)
  value <- value * effect


  # Apply missingness
  if (runif(1) < miss_prob) {
    value <- NA
  }

  return(value)
}
```

**Effect application rules:**
- `effect_size` applies only to `drug` treatment at `24h` timepoint
- All other conditions: effect = 1.0 (no change)

### 1.4 Expected Data Characteristics

| Metric | Expected Value |
|--------|----------------|
| Total observations | 500 × 24 = 12,000 |
| Missing rate | ~10-12% overall |
| MNAR correlation | r ≈ -0.4 (abundance vs missingness) |
| True positives | 50 peptides (drug effect at 24h) |
| Effect sizes | 1.8-2.5 fold |

### 1.5 Output Format

CSV with columns:
- `peptide`: Peptide identifier
- `gene_id`: Gene identifier (5 peptides per gene)
- `treatment`: "ctrl" or "drug"
- `timepoint`: "0h", "6h", or "24h"
- `bio_rep`: 1-4
- `abundance`: Numeric value (or NA)

---

## 2. Analysis Workflow

### 2.1 Data Import

```r
library(pepdiff)

dat <- read_pepdiff(
  "simulated_phospho.csv",
  id = "peptide",
  gene = "gene_id",
  value = "abundance",
  factors = c("treatment", "timepoint"),
  replicate = "bio_rep"
)

print(dat)
# Shows: 500 peptides, 24 samples, ~11% missing
```

### 2.2 Primary Analysis: GLM with Stratification

```r
results_glm <- compare(
  dat,
  compare = "treatment",
  ref = "ctrl",
  within = "timepoint",
  method = "glm"
)

print(results_glm)
# Shows: 3 comparisons (one per timepoint), n significant at each
```

**Expected results:**
- 0h: ~0-2 significant (false positives only)
- 6h: ~0-2 significant (false positives only)
- 24h: ~40-45 significant (true positives + few false positives)

### 2.3 Comparison: Conventional Workflow (Perseus-style)

Many proteomics tools use a conventional workflow:
1. Log2-transform abundances
2. Impute missing values (downshifted normal distribution)
3. Apply t-test per comparison

We implement this workflow using standard R functions to enable direct comparison:

```r
# ============================================================
# CONVENTIONAL WORKFLOW (mimics Perseus/standard approach)
# ============================================================

# Step 1: Log2 transform
log_data <- dat$data %>%
  mutate(log2_abundance = log2(value))

# Step 2: Impute missing values (downshifted normal)
# Perseus default: mean - 1.8*sd, width 0.3*sd
impute_downshifted <- function(x) {
  observed <- x[!is.na(x)]
  if (length(observed) < 2) return(x)

  mu_impute <- mean(observed) - 1.8 * sd(observed)
  sd_impute <- 0.3 * sd(observed)

  x[is.na(x)] <- rnorm(sum(is.na(x)), mean = mu_impute, sd = sd_impute)
  return(x)
}

imputed_data <- log_data %>%
  group_by(peptide) %>%
  mutate(log2_imputed = impute_downshifted(log2_abundance)) %>%
  ungroup()

# Step 3: t-test per peptide at 24h timepoint
conventional_results <- imputed_data %>%
  filter(timepoint == "24h") %>%
  group_by(peptide) %>%
  summarise(
    mean_ctrl = mean(log2_imputed[treatment == "ctrl"]),
    mean_drug = mean(log2_imputed[treatment == "drug"]),
    log2_fc = mean_drug - mean_ctrl,
    p_value = t.test(
      log2_imputed[treatment == "drug"],
      log2_imputed[treatment == "ctrl"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    significant = fdr < 0.05
  )
```

**Problems with conventional workflow:**

1. **Imputation bias with MNAR data**: Downshifted imputation assumes missingness is random within the low-abundance range. With true MNAR (missing = below detection), imputed values may be too high, inflating abundance estimates for low-abundance peptides.

2. **Normality assumption after log-transform**: t-test assumes normality. Log-transformed Gamma data is approximately normal, but heavy tails remain for some peptides.

3. **Variance homogeneity**: t-test assumes equal variance. Proteomics data typically shows mean-variance relationship even after log-transform.

4. **Information loss from imputation**: Imputed values treated as real data, ignoring uncertainty.

### 2.4 Comparison: Complete Cases Only

An alternative conservative approach - analyse only peptides with complete data:

```r
# Complete cases analysis (no imputation)
complete_peptides <- dat$data %>%
  filter(timepoint == "24h") %>%
  group_by(peptide) %>%
  summarise(n_complete = sum(!is.na(value))) %>%
  filter(n_complete == 8) %>%  # 4 ctrl + 4 drug

  pull(peptide)

complete_results <- dat$data %>%
  filter(peptide %in% complete_peptides, timepoint == "24h") %>%
  group_by(peptide) %>%
  summarise(
    log2_fc = log2(mean(value[treatment == "drug"]) /
                   mean(value[treatment == "ctrl"])),
    p_value = wilcox.test(
      value[treatment == "drug"],
      value[treatment == "ctrl"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    significant = fdr < 0.05
  )
```

**Problem**: Excludes low-abundance peptides (the ones with most missingness), biasing results toward high-abundance peptides and missing genuine biology.

### 2.5 Comparison Analysis: Pairwise Wilcoxon (pooled)

```r
# Unstratified pairwise comparison (ignoring timepoint structure)
results_pairwise <- compare(
  dat,
  compare = "treatment",
  ref = "ctrl",
  method = "pairwise",
  test = "wilcoxon"
)
```

**Expected results:**
- Fewer true positives detected (signal diluted across timepoints)
- Less precise effect estimates

### 2.7 Diagnostics

```r
diag <- plot_fit_diagnostics(results_glm)

# Check flagged peptides
diag$flagged
# Should show ~5-15 peptides with poor GLM fit

# View diagnostic plots
diag$plot
```

### 2.8 ART Comparison for Problematic Peptides

```r
# Rerun with ART for comparison
results_art <- compare(
  dat,
  compare = "treatment",
  ref = "ctrl",
  within = "timepoint",
  method = "art"
)
```

---

## 3. Figure 1 Panels

### 3.1 Panel A: Workflow Diagram

**Type**: Conceptual diagram (not data-driven)

**Content**:
```
┌─────────┐    ┌───────────────┐    ┌─────────────┐    ┌─────────────────┐
│   CSV   │ → │ read_pepdiff()│ → │ pepdiff_data│ → │    compare()    │
└─────────┘    └───────────────┘    └─────────────┘    └─────────────────┘
                                                               │
                                                               ▼
┌─────────┐    ┌───────────────┐    ┌─────────────────────────────────────┐
│  plots  │ ← │    plot()     │ ← │          pepdiff_results            │
└─────────┘    └───────────────┘    └─────────────────────────────────────┘
```

**Style**:
- Use pepdiff color scheme (cream #FFFFCC, orange #FD8D3C, red #BD0026)
- Clean boxes with rounded corners
- Arrows showing data flow

**Creation**: Manual in Inkscape or via DiagrammeR

### 3.2 Panel B: Volcano Plot

**Type**: Data-driven (ggplot2)

**Data source**: `results_glm` filtered to 24h timepoint

**Code**:
```r
library(ggplot2)

# Extract 24h results
res_24h <- results_glm$results %>%
  filter(comparison == "drug_vs_ctrl_24h")

# Add ground truth for coloring
res_24h <- res_24h %>%
  left_join(peptide_params %>% select(peptide, true_effect))

panel_b <- ggplot(res_24h, aes(x = log2_fc, y = -log10(fdr))) +
  geom_point(aes(color = significant & true_effect), alpha = 0.6) +
  scale_color_manual(
    values = c("grey70", "#BD0026"),
    labels = c("Not significant / FP", "True positive")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  labs(
    x = "log2 Fold Change (drug vs ctrl)",
    y = "-log10 FDR",
    title = "Treatment effect at 24h"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
```

**Key features**:
- Clear cluster of true positives (red) in upper right
- Few false positives
- Reference lines at FDR = 0.05 and |log2FC| = 1

### 3.3 Panel C: Method Comparison

**Type**: Data-driven (ggplot2)

**Data source**: Calculated metrics from all four analysis approaches

**Code**:
```r
# Calculate metrics for each method
calc_metrics <- function(results_df, ground_truth, method_name) {
  # Join with ground truth
  merged <- results_df %>%
    left_join(ground_truth %>% select(peptide, true_effect), by = "peptide")

  tibble(
    method = method_name,
    true_positives = sum(merged$significant & merged$true_effect, na.rm = TRUE),
    false_positives = sum(merged$significant & !merged$true_effect, na.rm = TRUE),
    false_negatives = sum(!merged$significant & merged$true_effect, na.rm = TRUE),
    sensitivity = true_positives / (true_positives + false_negatives),
    fdr = false_positives / (true_positives + false_positives)
  )
}

# Gather metrics from all methods
comparison_metrics <- bind_rows(
 calc_metrics(results_glm$results %>% filter(comparison == "drug_vs_ctrl_24h"),
               peptide_params, "pepdiff GLM"),
  calc_metrics(conventional_results, peptide_params, "Conventional"),
  calc_metrics(complete_results, peptide_params, "Complete cases"),
  calc_metrics(results_pairwise$results, peptide_params, "Pairwise pooled")
)

# Reshape for plotting
plot_data <- comparison_metrics %>%
  select(method, sensitivity, fdr) %>%
  pivot_longer(cols = c(sensitivity, fdr), names_to = "metric", values_to = "value") %>%
  mutate(
    method = factor(method, levels = c("pepdiff GLM", "Conventional",
                                        "Complete cases", "Pairwise pooled")),
    metric = factor(metric, levels = c("sensitivity", "fdr"),
                    labels = c("Sensitivity", "FDR"))
  )

panel_c <- ggplot(plot_data, aes(x = method, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Sensitivity" = "#FD8D3C", "FDR" = "#BD0026")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "Value",
    fill = NULL,
    title = "Method comparison"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
```

**Key features**:
- pepdiff GLM shows highest sensitivity (~86%) with lowest FDR (~5%)
- Conventional approach has inflated FDR (~15%) due to imputation bias
- Complete cases is conservative (low FDR) but loses sensitivity (~54%)
- Pairwise pooled dilutes signal, losing sensitivity

**Message**: pepdiff achieves the best balance of sensitivity and specificity.

### 3.4 Panel D: Diagnostic QQ Plots

**Type**: Data-driven (ggplot2)

**Data source**: `diag$flagged` and residuals from `results_glm$diagnostics`

**Code**:
```r
# Get one good-fit and one poor-fit peptide
good_pep <- results_glm$diagnostics %>%
  filter(converged, deviance < 1) %>%
  slice(1) %>%
  pull(peptide)

poor_pep <- diag$flagged %>%
  slice(1) %>%
  pull(peptide)

# Extract standardized residuals
good_resid <- results_glm$diagnostics %>%
  filter(peptide == good_pep) %>%
  pull(std_residuals) %>%
  unlist()

poor_resid <- results_glm$diagnostics %>%
  filter(peptide == poor_pep) %>%
  pull(std_residuals) %>%
  unlist()

# Create QQ data
qq_data <- bind_rows(
  tibble(residual = good_resid, type = "Good fit"),
  tibble(residual = poor_resid, type = "Poor fit")
)

panel_d <- ggplot(qq_data, aes(sample = residual)) +
  stat_qq() +
  stat_qq_line(color = "#BD0026") +
  facet_wrap(~type) +
  labs(
    x = "Theoretical quantiles",
    y = "Sample quantiles",
    title = "Model diagnostics"
  ) +
  theme_minimal()
```

**Key features**:
- Left panel: Points follow diagonal (Gamma GLM appropriate)
- Right panel: S-shaped deviation (heavy tails, consider ART)
- Clear visual difference guides method choice

### 3.5 Supplementary Figure: Heatmap of Temporal Dynamics

**Type**: Data-driven (ComplexHeatmap or ggplot2)

**Data source**: Significant peptides from `results_glm`, abundance matrix

**Code**:
```r
# Get significant peptides at 24h
sig_peptides <- results_glm$results %>%
  filter(significant, comparison == "drug_vs_ctrl_24h") %>%
  pull(peptide)

# Create abundance matrix
abundance_matrix <- dat$data %>%
  filter(peptide %in% sig_peptides) %>%
  group_by(peptide, treatment, timepoint) %>%
  summarise(mean_abundance = mean(value, na.rm = TRUE), .groups = "drop") %>%
  unite("condition", treatment, timepoint) %>%
  pivot_wider(names_from = condition, values_from = mean_abundance) %>%
  column_to_rownames("peptide") %>%
  as.matrix()

# Z-score normalize rows
abundance_z <- t(scale(t(abundance_matrix)))

# Order columns logically
col_order <- c("ctrl_0h", "ctrl_6h", "ctrl_24h", "drug_0h", "drug_6h", "drug_24h")
abundance_z <- abundance_z[, col_order]

# Heatmap with ggplot2
heatmap_data <- abundance_z %>%
  as.data.frame() %>%
  rownames_to_column("peptide") %>%
  pivot_longer(-peptide, names_to = "condition", values_to = "z_score")

heatmap_data$condition <- factor(heatmap_data$condition, levels = col_order)

supp_heatmap <- ggplot(heatmap_data, aes(x = condition, y = peptide, fill = z_score)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#BD0026",
    midpoint = 0, name = "Z-score"
  ) +
  labs(
    x = "Condition",
    y = "Peptide",
    title = "Significant peptides across timepoints"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
```

**Key features**:
- Rows: ~40-45 significant peptides (clustered)
- Columns: 6 conditions in logical order
- Clear pattern: effect emerges at drug_24h
- Color: Blue (low) → White → Red (high)

---

## 4. Validation Metrics

### 4.1 Primary Analysis Performance (pepdiff GLM)

| Metric | Expected | Acceptable Range |
|--------|----------|------------------|
| True positives (24h) | 42-45 | 38-48 |
| False positives (24h) | 2-3 | 0-5 |
| False positives (0h + 6h) | 2-4 | 0-6 |
| FDR (24h, empirical) | 0.05-0.07 | < 0.10 |
| Sensitivity | 0.84-0.90 | > 0.75 |

### 4.2 Method Comparison Summary

| Metric | pepdiff GLM | Conventional (impute + t-test) | Complete cases | Pairwise pooled |
|--------|-------------|-------------------------------|----------------|-----------------|
| Peptides tested | 500 | 500 | ~350-400 | 500 |
| True positives | 42-45 | 35-40 | 25-30 | 25-32 |
| False positives | 2-3 | 5-10 | 1-2 | 3-5 |
| Sensitivity | ~0.86 | ~0.74 | ~0.55 | ~0.56 |
| FDR (empirical) | ~0.05 | ~0.15 | ~0.05 | ~0.12 |

**Key insights:**

1. **pepdiff GLM vs Conventional**: pepdiff detects ~15% more true positives with lower FDR. Imputation inflates false positives because imputed values for low-abundance peptides (which have MNAR missingness) are too high, creating artificial "effects."

2. **pepdiff GLM vs Complete cases**: Complete cases analysis is conservative (low FDR) but misses ~35% of true positives because low-abundance peptides with genuine effects are excluded due to missingness.

3. **pepdiff GLM vs Pairwise pooled**: Stratified analysis captures the timepoint-specific effect; pooled analysis dilutes the signal across all timepoints, losing sensitivity.

### 4.3 Conventional Workflow Problems (Detailed)

| Issue | Consequence | pepdiff Solution |
|-------|-------------|------------------|
| Imputation with MNAR | False positives from inflated low-abundance values | No imputation; Gamma GLM handles incomplete cases |
| Log-transform + t-test | Residual non-normality; variance heterogeneity | Native Gamma distribution; mean-variance relationship modeled |
| Separate per-timepoint tests | Correct approach but no framework for factorial structure | `within` parameter for stratified contrasts from single model |
| No diagnostics | Can't identify problematic peptides | `plot_fit_diagnostics()` flags poor fits |

### 4.4 Diagnostic Flagging

| Metric | Expected |
|--------|----------|
| Peptides flagged | 5-15 |
| Flagged with true effect | 1-3 |
| ART recovers flagged TPs | Yes |

---

## 5. Supplementary Code Structure

### 5.1 File: `supplementary_example.Rmd`

```yaml
---
title: "pepdiff: Supplementary Example Code"
output:
  pdf_document:
    toc: true
---
```

### 5.2 Sections

1. **Setup**: Load packages, set seed
2. **Data Generation**: Create simulated dataset with ground truth
3. **Data Import**: `read_pepdiff()` demonstration
4. **Primary Analysis**: GLM with stratification
5. **Conventional Workflow**: Imputation + t-test (Perseus-style)
6. **Complete Cases Analysis**: Conservative alternative
7. **Comparison: GLM vs Pairwise**: Side-by-side performance
8. **Method Comparison Table**: Summary of all approaches
9. **Diagnostics**: `plot_fit_diagnostics()` demonstration
10. **Comparison: GLM vs ART**: For problematic peptides
11. **Figure 1**: Generate all panels
12. **Figure 2** (Supplementary): Method comparison bar chart
13. **Session Info**: Reproducibility

### 5.3 Runtime Expectations

| Step | Expected Time |
|------|---------------|
| Data generation | < 1 sec |
| Data import | < 1 sec |
| GLM analysis (500 peptides) | 30-60 sec |
| ART analysis (500 peptides) | 45-90 sec |
| Diagnostics | < 5 sec |
| Figure generation | < 10 sec |
| **Total** | **~2-3 min** |

---

## 6. Key Code Snippets for Paper

### Box 1: Basic Workflow

```r
library(pepdiff)

# Import data
dat <- read_pepdiff(
  "experiment.csv",
  id = "peptide",
  gene = "gene_id",
  value = "abundance",
  factors = c("treatment", "timepoint"),
  replicate = "bio_rep"
)

# Analyze with stratification
results <- compare(
  dat,
  compare = "treatment",
  ref = "ctrl",
  within = "timepoint",
  method = "glm"
)

# Visualize
plot(results)
```

### Box 2: Diagnostic Check

```r
# Assess model fit
diag <- plot_fit_diagnostics(results)

# View flagged peptides
diag$flagged

# If many flagged, consider ART
results_art <- compare(dat, compare = "treatment",
                       ref = "ctrl", method = "art")
```

---

## 7. Dependencies

```r
# Required packages
library(pepdiff)      # This package
library(tidyverse)    # Data manipulation
library(cowplot)      # Figure composition
library(ComplexHeatmap)  # Optional: better heatmaps

# For paper rendering
library(rmarkdown)
library(knitr)
```

---

## 8. Output Files

| File | Description |
|------|-------------|
| `simulated_phospho.csv` | Generated dataset |
| `peptide_params.csv` | Ground truth parameters |
| `figure_1.pdf` | Combined 4-panel figure |
| `panel_a_workflow.svg` | Workflow diagram (vector) |
| `supplementary_example.pdf` | Rendered supplementary |
