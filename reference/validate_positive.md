# Validate positive values for Gamma GLM

Checks that all non-NA values are strictly positive (\> 0). Gamma GLM
requires positive values.

## Usage

``` r
validate_positive(values, name = "values")
```

## Arguments

- values:

  Numeric vector of values to check

- name:

  Name of the variable for error messages

## Value

TRUE invisibly if valid, otherwise throws an error
