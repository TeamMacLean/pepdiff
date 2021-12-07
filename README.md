
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pepdiff

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/TeamMacLean/pepdiff.svg?branch=master)](https://travis-ci.org/TeamMacLean/pepdiff)
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
