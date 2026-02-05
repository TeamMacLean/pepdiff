# Run statistical model for all peptides

Applies GLM or ART model to each peptide in the dataset.

## Usage

``` r
run_models(data, compare, ref = NULL, method = "glm")
```

## Arguments

- data:

  A pepdiff_data object

- compare:

  Factor to compare

- ref:

  Reference level

- method:

  One of "glm" or "art"

## Value

A list with results for each peptide
