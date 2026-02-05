# Build results tibble from pairwise tests

Build results tibble from pairwise tests

## Usage

``` r
build_pairwise_results(
  test_results,
  data,
  compare,
  ref,
  test,
  alpha,
  fdr_method
)
```

## Arguments

- test_results:

  List of test results for each peptide

- data:

  Original pepdiff_data object

- compare:

  Factor being compared

- ref:

  Reference level

- test:

  Test method name

- alpha:

  Significance threshold

- fdr_method:

  FDR correction method

## Value

A tibble in long format
