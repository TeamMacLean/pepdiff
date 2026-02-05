
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pepdiff <a href="https://teammaclean.github.io/pepdiff/"><img src="man/figures/logo.png" align="right" height="138" alt="pepdiff website" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/pepdiff)](https://CRAN.R-project.org/package=pepdiff)
[![R-CMD-check](https://github.com/TeamMacLean/pepdiff/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TeamMacLean/pepdiff/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/TeamMacLean/pepdiff/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/TeamMacLean/pepdiff/actions/workflows/pkgdown.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub
issues](https://img.shields.io/github/issues/TeamMacLean/pepdiff)](https://github.com/TeamMacLean/pepdiff/issues)
[![GitHub
stars](https://img.shields.io/github/stars/TeamMacLean/pepdiff)](https://github.com/TeamMacLean/pepdiff/stargazers)
<!-- badges: end -->

**Differential abundance analysis for PRM proteomics data.**

pepdiff helps proteomics researchers answer: “What’s differentially
abundant?”

## Features

- **GLM analysis** – Gamma GLM with emmeans for factorial designs
- **ART analysis** – Non-parametric alternative for heavy-tailed data
- **Pairwise tests** – Wilcoxon, bootstrap-t, Bayes factor, rank
  products
- **Stratified comparisons** – Analyse effects within factor levels
- **Fit diagnostics** – Visual checks for GLM model assumptions
- **Rich visualizations** – Volcano plots, heatmaps, PCA, p-value
  histograms

## Installation

Install from GitHub:

``` r
# install.packages("pak")
pak::pak("TeamMacLean/pepdiff")
```

## Quick Start

``` r
library(pepdiff)

# Import data
dat <- read_pepdiff(
  "data.csv",
  id = "peptide",
  gene = "gene_id",
  value = "abundance",
  factors = c("treatment", "timepoint"),
  replicate = "bio_rep"
)

# Run differential analysis
results <- compare(
  dat,
  compare = "treatment",
  ref = "ctrl",
  method = "glm"
)

# Visualize results
plot(results)
```

## Documentation

- **[Getting
  Started](https://teammaclean.github.io/pepdiff/articles/basic_workflow.html)**
  – Basic workflow
- **[GLM
  Analysis](https://teammaclean.github.io/pepdiff/articles/glm_analysis.html)**
  – Factorial designs with GLM
- **[ART
  Analysis](https://teammaclean.github.io/pepdiff/articles/art_analysis.html)**
  – Non-parametric alternative
- **[Checking Model
  Fit](https://teammaclean.github.io/pepdiff/articles/checking_fit.html)**
  – Diagnostic plots
- **[Pairwise
  Tests](https://teammaclean.github.io/pepdiff/articles/pairwise_tests.html)**
  – Direct two-group comparisons
- **[Function
  Reference](https://teammaclean.github.io/pepdiff/reference/index.html)**
  – Full API

## Companion Package

**peppwR** answers “How many samples do I need?” (power analysis)
**pepdiff** answers “What’s differentially abundant?” (analysis)

See [peppwR](https://github.com/TeamMacLean/peppwR) for experimental
design planning.

## Workflow

``` mermaid
flowchart LR
    A[CSV] --> B[read_pepdiff]
    B --> C[pepdiff_data]
    C --> D[compare]
    D --> E[pepdiff_results]
    E --> F[plot]

    style A fill:#FFFFCC,stroke:#BD0026
    style B fill:#FD8D3C,stroke:#BD0026,color:#fff
    style C fill:#FFFFCC,stroke:#BD0026
    style D fill:#FD8D3C,stroke:#BD0026,color:#fff
    style E fill:#FFFFCC,stroke:#BD0026
    style F fill:#FD8D3C,stroke:#BD0026,color:#fff
```

## Citation

If you use pepdiff in your research, please cite:

    MacLean, D. (2026). pepdiff: Differential Abundance Analysis for PRM
    Proteomics Data. R package version 1.0.0.
    https://github.com/TeamMacLean/pepdiff

## Contributing

Contributions welcome! Please open an
[issue](https://github.com/TeamMacLean/pepdiff/issues) or submit a pull
request.

## License

MIT
