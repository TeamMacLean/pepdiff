
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pepdiff

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/TeamMacLean/pepdiff/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TeamMacLean/pepdiff/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of pepdiff is to provide functions and plots for normalising
peptide data from ‘PRM’ experiments.

## Installation

Install direct from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("TeamMacLean/pepdiff")
```

Or via `renv`.

``` r
# install.packages("devtools")
renv::install("TeamMacLean/pepdiff)
```

## Use

The package contains vignettes describing the basic workflows and you
should read them before using the package.

There are four vignettes and they should be read in order as they follow
the analysis from start to finish.

### Data Checks

How to check the coverage in the data and assess things like it’s
distribution.

``` r
library(pepdiff)
vignette('data-checks', package="pepdiff")
```

### Doing Comparisons

How to perform comparisons with selected tests

``` r
library(pepdiff)
vignette('doing-comparisons', package="pepdiff")
```

### Check Results

How to evaluate the tests you’ve chosen and make wiser interpretations
and experimental designs

``` r
library(pepdiff)
vignette('check-results', package="pepdiff")
```

### Analysis

How to perform an analysis of a broad set of comparisons

``` r
library(pepdiff)
vignette('analysis', package="pepdiff")
```
