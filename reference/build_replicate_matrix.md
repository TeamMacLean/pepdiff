# Build replicate matrix from pepdiff_data

Pivots data to matrix format (peptides × replicates) for a given factor
level.

## Usage

``` r
build_replicate_matrix(data, compare, level)
```

## Arguments

- data:

  pepdiff_data object

- compare:

  Factor column name

- level:

  Level of the factor to extract

## Value

Matrix with peptides as rows and replicates as columns
