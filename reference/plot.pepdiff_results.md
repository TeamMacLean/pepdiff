# Plot method for pepdiff_results

Creates a multi-panel plot showing volcano, p-value/BF histogram, and FC
distribution. Automatically dispatches to BF-specific plots when results
are from bayes_t test.

## Usage

``` r
# S3 method for class 'pepdiff_results'
plot(x, ...)
```

## Arguments

- x:

  A pepdiff_results object

- ...:

  Additional arguments (ignored)

## Value

A cowplot grid of plots
