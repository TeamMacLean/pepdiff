# Compare within strata

Runs comparisons separately within each level of stratifying factor(s).

## Usage

``` r
compare_within_strata(
  data,
  compare,
  ref,
  within,
  method,
  test,
  alpha,
  fdr_method,
  bf_threshold = 3
)
```

## Arguments

- data:

  pepdiff_data object

- compare:

  Factor to compare

- ref:

  Reference level

- within:

  Factor(s) to stratify by

- method:

  Analysis method

- test:

  Pairwise test (if applicable)

- alpha:

  Significance threshold

- fdr_method:

  FDR correction method

## Value

Combined results across strata
