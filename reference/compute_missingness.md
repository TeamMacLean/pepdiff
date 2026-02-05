# Compute missingness statistics per peptide

Calculates the NA rate and MNAR (Missing Not At Random) score for each
peptide. MNAR score is based on the relationship between missingness and
mean abundance.

## Usage

``` r
compute_missingness(data)
```

## Arguments

- data:

  A data frame with 'peptide' and 'value' columns

## Value

A tibble with columns: peptide, na_rate, mnar_score, mean_abundance
