# Validate factor columns exist in data

Checks that all specified factor columns exist in the data frame and are
not empty.

## Usage

``` r
validate_factors(data, factors)
```

## Arguments

- data:

  A data frame

- factors:

  Character vector of factor column names

## Value

TRUE invisibly if valid, otherwise throws an error
