# Bayes factor distribution plot

Creates a histogram of log10(BF) values with reference lines at standard
thresholds.

## Usage

``` r
plot_bf_distribution(results, comparison = NULL)
```

## Arguments

- results:

  A pepdiff_results object from bayes_t test

- comparison:

  Optional comparison to filter by

## Value

A ggplot object
