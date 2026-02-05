# Apply FDR correction within groups

Applies Benjamini-Hochberg FDR correction within each comparison group.

## Usage

``` r
apply_fdr_by_comparison(results, method = "BH")
```

## Arguments

- results:

  A data frame with 'comparison' and 'p_value' columns

- method:

  FDR correction method (default "BH")

## Value

Data frame with added 'fdr' column
