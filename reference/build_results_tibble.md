# Build results tibble from model output

Combines model results into the standard long-format results tibble.

## Usage

``` r
build_results_tibble(model_results, data, compare, method, alpha, fdr_method)
```

## Arguments

- model_results:

  List of results from run_models()

- data:

  Original pepdiff_data object

- compare:

  Factor being compared

- method:

  Method used

- alpha:

  Significance threshold

- fdr_method:

  FDR correction method

## Value

A tibble in long format
