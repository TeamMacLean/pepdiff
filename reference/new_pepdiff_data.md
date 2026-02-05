# Create a new pepdiff_data object

Low-level constructor for pepdiff_data objects. Use \[read_pepdiff()\]
for user-facing data import.

## Usage

``` r
new_pepdiff_data(data, factors, design, missingness, peptides, call)
```

## Arguments

- data:

  A tibble with peptide, gene_id, factor columns, bio_rep, value

- factors:

  Character vector of factor column names

- design:

  Tibble of unique factor combinations with counts

- missingness:

  Tibble of peptide missingness statistics

- peptides:

  Character vector of unique peptide IDs

- call:

  The original function call

## Value

A pepdiff_data object
