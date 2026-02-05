# Compare using Rank Products (full matrix approach)

Rank Products requires ranking across ALL peptides simultaneously, not
within single peptides. This function handles the matrix-based approach
using the RankProd package.

## Usage

``` r
compare_pairwise_rankprod(data, compare, ref, alpha, fdr_method)
```

## Arguments

- data:

  pepdiff_data object

- compare:

  Factor to compare

- ref:

  Reference level

- alpha:

  Significance threshold

- fdr_method:

  FDR correction method

## Value

List with results and diagnostics tibbles
