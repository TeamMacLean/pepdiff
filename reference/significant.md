# Extract significant results

Extract significant results

## Usage

``` r
significant(x, alpha = NULL, by_fdr = TRUE, bf_threshold = NULL)
```

## Arguments

- x:

  A pepdiff_results object

- alpha:

  Significance threshold for p-value tests (default uses analysis alpha)

- by_fdr:

  Logical, use FDR-adjusted p-values for p-value tests (default TRUE)

- bf_threshold:

  BF threshold for bayes_t results (default uses analysis bf_threshold)

## Value

A tibble of significant results
