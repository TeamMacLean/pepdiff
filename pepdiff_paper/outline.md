# pepdiff: Differential abundance analysis for phosphoproteomics experiments

## Target Journal
**Bioinformatics** (Oxford) - Application Note

## Format Requirements
- ~2,000 words + 1 figure (max 4 pages / 2,600 words)
- Software name in title
- Supplementary material allowed

## Companion Framing

**peppwR** answers: "How many samples do I need?" (experimental design)
**pepdiff** answers: "What's differentially abundant?" (analysis)

Together they provide an end-to-end workflow for phosphoproteomics experiments.

---

## Abstract (~150 words)

**Background:** Phosphoproteomics experiments generate complex factorial designs requiring appropriate statistical models for differential abundance analysis. Existing tools either assume simple two-group comparisons or require extensive statistical expertise to specify models correctly.

**Results:** We present pepdiff, an R package for differential abundance analysis of phosphoproteomics data. pepdiff implements Gamma GLM with emmeans-based contrasts for factorial designs, Aligned Rank Transform (ART) for non-parametric alternatives, and four pairwise tests (Wilcoxon, bootstrap-t, Bayes factor, rank products) for simple comparisons. Built-in diagnostics help users assess model fit and choose appropriate methods. A simple interface handles common designs while a formula interface enables complex contrasts. Stratified comparisons allow analysis within factor levels.

**Availability:** pepdiff is freely available from GitHub (https://github.com/TeamMacLean/pepdiff). Documentation at https://teammaclean.github.io/pepdiff/. Companion package peppwR provides power analysis for experimental design.

---

## 1. Introduction (~400 words)

### The analysis gap in proteomics

- Phosphoproteomics experiments increasingly use factorial designs (treatment × timepoint, genotype × condition)
- Simple pairwise tests ignore factorial structure, losing power and missing interactions
- Proper GLM-based analysis requires statistical expertise to specify correctly
- Gap between simple tools (inflexible) and complex tools (inaccessible)

### Why proteomics data needs special handling

- **Right-skewed distributions**: Abundance data better modeled by Gamma than Normal
- **Heterogeneous variance**: Variance scales with mean (heteroscedasticity)
- **Missing data**: Systematic missingness at low abundances
- **Multiple testing**: Thousands of peptides require FDR correction

### Existing tools and their limitations

| Tool | Limitation |
|------|------------|
| Perseus | Powerful but complex; steep learning curve for factorial designs |
| limma | Designed for microarrays; assumes log-normal after transformation |
| MSstats | Comprehensive but complex; focuses on protein-level inference |
| Simple t-tests | Ignore factorial structure; assume normality |

### pepdiff addresses these gaps

- Gamma GLM naturally handles right-skewed, heteroscedastic data
- Simple interface for common designs, formula interface for complex ones
- Built-in diagnostics guide method choice (GLM vs ART)
- Stratified comparisons for within-factor analysis
- Companion to peppwR for complete experimental workflow

---

## 2. Implementation (~600 words)

### 2.1 Workflow overview

```
CSV → read_pepdiff() → pepdiff_data → compare() → pepdiff_results → plots
```

S3 classes provide consistent interface:
- `pepdiff_data`: Import, validation, missingness assessment
- `pepdiff_results`: Results, diagnostics, visualization methods

### 2.2 Three analysis methods

| Method | Model | Use Case |
|--------|-------|----------|
| `"glm"` (default) | Gamma GLM + emmeans | Factorial designs, skewed data |
| `"art"` | Aligned Rank Transform | Non-parametric alternative |
| `"pairwise"` | Direct two-group tests | Simple comparisons |

### 2.3 Gamma GLM approach

- Gamma distribution: Natural for positive, right-skewed abundance data
- Log link: Multiplicative effects on original scale
- emmeans: Extract pairwise contrasts from fitted model
- Per-peptide fitting: Each peptide modeled independently

**Why Gamma, not log-transform + linear model?**
- Handles zeros and near-zeros without arbitrary constants
- Variance structure (mean-variance relationship) modeled correctly
- Inference on original scale, not transformed

### 2.4 Aligned Rank Transform (ART)

- Non-parametric alternative when GLM assumptions violated
- Preserves factorial structure (unlike Kruskal-Wallis)
- Uses ARTool package [@wobbrock2011art]
- Same emmeans-based contrast extraction

### 2.5 Pairwise tests

Four tests matching peppwR for workflow consistency:
- **Wilcoxon rank-sum**: Non-parametric, robust
- **Bootstrap-t**: Handles non-normality via resampling
- **Bayes factor t-test**: Bayesian evidence quantification
- **Rank products**: Non-parametric, handles heterogeneity [@breitling2004rankprod]

### 2.6 Model diagnostics

`plot_fit_diagnostics()` provides 4-panel assessment:
1. **Deviance distribution**: Overall model fit quality
2. **Deviance vs fold change**: Detect effect-size-dependent misfit
3. **Sample residual plots**: Individual peptide diagnostics
4. **Pooled QQ plot**: Distributional assumption check

Flagged peptides: High deviance suggests poor fit; consider ART method.

### 2.7 Stratified comparisons

`within` parameter enables comparisons within factor levels:
```r
compare(data, compare = "treatment", ref = "ctrl", within = "timepoint")
```
Produces separate treatment effects at each timepoint.

### 2.8 Simple and formula interfaces

**Simple interface** (most users):
```r
compare(data, compare = "treatment", ref = "ctrl", method = "glm")
```

**Formula interface** (power users):
```r
compare(data, contrast = pairwise ~ treatment | timepoint, method = "glm")
```

---

## 3. Example Application (~450 words)

### 3.1 Simulated phosphoproteomics dataset (~80 words)

- Generate realistic factorial design: treatment (ctrl, drug) × timepoint (0h, 6h, 24h)
- 500 peptides with heterogeneous Gamma-distributed abundances
- 4 biological replicates per condition (24 samples total)
- ~10% missingness with MNAR pattern
- 50 peptides with true treatment effect (2-fold at 24h)

### 3.2 GLM analysis of factorial design (~100 words)

```r
dat <- read_pepdiff("experiment.csv", ...)
results <- compare(dat, compare = "treatment", ref = "ctrl",
                   within = "timepoint", method = "glm")
```

**Primary result**: "42 of 50 true positives detected at FDR < 0.05; 3 false positives"

Key insight: Effect emerges at 24h timepoint; stratified analysis reveals temporal dynamics that pooled analysis would miss.

### 3.3 Comparison 1: pepdiff vs conventional workflow (~120 words)

**Key comparison demonstrating why pepdiff's approach matters:**

Conventional proteomics workflow (as used in Perseus and similar tools):
1. Log2-transform data
2. Impute missing values (downshifted normal)
3. Apply t-test

Results comparison at 24h timepoint:

| Approach | True Positives | False Positives | Sensitivity |
|----------|---------------|-----------------|-------------|
| pepdiff GLM | 43 | 2 | 0.86 |
| Conventional | 37 | 8 | 0.74 |
| Complete cases only | 27 | 1 | 0.54 |

**Message**: Imputation with MNAR data inflates false positives. pepdiff's Gamma GLM handles missing data naturally without imputation, improving both sensitivity and specificity.

### 3.4 Comparison 2: GLM vs pairwise tests (~80 words)

- Pairwise Wilcoxon (ignoring timepoint): 28 true positives, higher FDR
- GLM with stratification: 43 true positives, controlled FDR

**Message**: Stratified analysis captures factorial structure; pooled analysis dilutes timepoint-specific effects.

### 3.5 Comparison 3: GLM vs ART for problematic peptides (~80 words)

Introduce subset with heavy-tailed distributions:
- GLM: Poor fit (high deviance), unreliable p-values
- ART: Robust results despite distributional violations

**Message**: Diagnostics identify peptides where GLM assumptions fail; ART provides a valid alternative.

### 3.6 Model diagnostics in action (~60 words)

`plot_fit_diagnostics()` reveals:
- Most peptides: QQ plot follows diagonal (Gamma appropriate)
- Flagged peptides: Systematic deviation suggests heavy tails
- Decision: Use GLM for majority, consider ART for flagged subset or entire dataset

### 3.7 Figure panels (~60 words)

- **Panel A**: Workflow diagram (data → compare → results → plots)
- **Panel B**: Volcano plot showing treatment effects at 24h timepoint
- **Panel C**: Method comparison bar chart (sensitivity/FDR by approach)
- **Panel D**: Diagnostic QQ plots (good vs poor GLM fit)

---

## 4. Discussion (~200 words)

### Key contributions

- First dedicated differential abundance tool for phosphoproteomics factorial designs
- Gamma GLM naturally handles proteomics data properties
- Built-in diagnostics guide method selection
- Stratified comparisons for temporal/factorial analysis
- Companion to peppwR for complete workflow

### Limitations

- Per-peptide models: No borrowing of information across peptides (unlike limma's empirical Bayes)
- Single-factor contrasts: Complex multi-way interactions require formula interface
- Computational cost: GLM fitting for thousands of peptides takes minutes

### Future directions

- Protein-level rollup with peptide-level inference
- Mixed models for repeated measures designs
- Integration with MSstats/Perseus workflows

---

## Figure 1: Multi-panel overview

**Panel A**: Workflow diagram (conceptual)
- Flow: `CSV → read_pepdiff() → compare() → plot()`
- Clean schematic showing S3 class structure

**Panel B**: Volcano plot (data-driven)
- X-axis: log2 fold change
- Y-axis: -log10 FDR
- Color: Significant (red) vs non-significant (grey)
- Shows treatment effect at 24h timepoint
- Message: Clear separation of true positives

**Panel C**: Method comparison (data-driven)
- Grouped bar chart or dot plot
- X-axis: Method (pepdiff GLM, Conventional, Complete cases, Pairwise)
- Y-axis: Metric value
- Two metrics: Sensitivity (filled) and FDR (hollow)
- Message: pepdiff achieves highest sensitivity with lowest FDR

**Panel D**: Diagnostic QQ plots (data-driven)
- Side-by-side QQ plots: Good fit vs poor fit peptide
- Left: Points follow diagonal (Gamma appropriate)
- Right: S-shaped deviation (heavy tails, consider ART)
- Message: Built-in diagnostics guide method choice

---

## Key Messages (3 main points)

1. **No imputation needed**: Gamma GLM handles missing data naturally. Conventional imputation (downshifted normal) inflates false positives with MNAR data—pepdiff improves sensitivity by ~15% while reducing FDR.

2. **Factorial structure matters**: Stratified comparisons capture experimental design; pooled pairwise tests dilute timepoint-specific effects and lose ~40% of true positives.

3. **Diagnostics guide decisions**: Built-in `plot_fit_diagnostics()` identifies peptides where GLM assumptions fail, guiding users to ART when appropriate.

---

## Supplementary Material

- Fully reproducible R code for all analyses and Figure 1
- Extended diagnostics examples
- Comparison with limma/MSstats on benchmark data
- Runtime benchmarks

---

## References (to collect)

- emmeans package [@lenth2022emmeans]
- ARTool package [@wobbrock2011art]
- Gamma GLM for count/abundance data [@ver2007quasi]
- Rank products [@breitling2004rankprod]
- Perseus [@tyanova2016perseus]
- MSstats [@choi2014msstats]
- limma [@ritchie2015limma]
- peppwR (companion package)
- Benjamini-Hochberg FDR [@benjamini1995fdr]

---

## File Structure

```
pepdiff_paper/
├── outline.md                   # This file
├── paper.Rmd                    # Main manuscript
├── references.bib               # BibTeX references
├── example_application_spec.md  # Spec for supplementary Rmd
├── supplementary_example.Rmd    # Reproducible code supplement
└── figures/
    ├── figure_1.pdf             # Combined figure for journal
    └── panel_a_workflow.svg     # Workflow diagram (vector)
```

---

## Rendering Pipeline

Same as peppwR:

```r
# Default PDF
rmarkdown::render("paper.Rmd")

# bioRxiv preprint
rmarkdown::render("paper.Rmd", output_format = rticles::arxiv_article())

# Bioinformatics journal
rmarkdown::render("paper.Rmd", output_format = rticles::bioinformatics_article())
```

---

## Author Contributions

D.M. conceived the project, designed the software architecture, implemented the package, and wrote the manuscript.

## Companion Package Statement

pepdiff is designed as a companion to peppwR [@maclean2026peppwr]. While peppwR addresses experimental design through power analysis ("How many samples do I need?"), pepdiff provides the subsequent differential abundance analysis ("What's differentially abundant?"). Together, they offer an end-to-end workflow for phosphoproteomics experiments.
