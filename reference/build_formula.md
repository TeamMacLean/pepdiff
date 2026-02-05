# Build formula from factor names

Creates a formula string for GLM fitting.

## Usage

``` r
build_formula(response, factors, interaction = TRUE)
```

## Arguments

- response:

  Name of response variable

- factors:

  Character vector of factor names

- interaction:

  Logical, include interactions? (default TRUE)

## Value

A formula object
