# Fit GLM and extract contrasts for one peptide

Convenience function that fits a GLM and extracts contrasts in one step.

## Usage

``` r
fit_and_extract_glm(data, response, factors, compare, ref, peptide_id)
```

## Arguments

- data:

  Data for one peptide

- response:

  Name of the response variable

- factors:

  Character vector of factor names

- compare:

  Which factor to compare

- ref:

  Reference level

- peptide_id:

  Peptide identifier

## Value

A list with model fit results and contrasts
