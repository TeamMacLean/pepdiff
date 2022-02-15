
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pepdiff

<!-- badges: start -->

[![R-CMD-check](https://github.com/TeamMacLean/pepdiff/workflows/R-CMD-check/badge.svg)](https://github.com/TeamMacLean/pepdiff/actions)
<!-- badges: end -->

The goal of pepdiff is to provide functions and plots for normalising
peptide data from ‘PRM’ experiments.

## Installation

Install direct from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("TeamMacLean/pepdiff")
```

Or via `renv` . Some bioconductor packages won’t automatically install
this way, so do them first.

``` r
renv::install("bioc::ComplexHeatmap)
renv::install("bioc::tidybulk)

renv::install("TeamMacLean/pepdiff)
```

## Use

Basic functions and a typical workflow are outlined in the `basic_use`
vignette.

``` r
library(pepdiff)
vignette('basic_use', package="pepdiff")
```
