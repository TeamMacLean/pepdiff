# Extract model diagnostics

Creates a summary tibble of model convergence and quality metrics.

## Usage

``` r
extract_diagnostics(model_results)
```

## Arguments

- model_results:

  List of model results from run_models()

## Value

A tibble with diagnostic information including residuals and fitted
values (list columns)
