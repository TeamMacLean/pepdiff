# Compute design summary from data

Creates a summary of factor combinations with replicate and peptide
counts.

## Usage

``` r
compute_design_summary(data, factors)
```

## Arguments

- data:

  A data frame with factor columns, 'bio_rep', and 'peptide' columns

- factors:

  Character vector of factor column names

## Value

A tibble with factor columns plus n_reps and n_peptides
