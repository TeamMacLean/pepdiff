# Fit a Gamma GLM for a single peptide

Fits a generalized linear model with Gamma family and log link to model
abundance values as a function of experimental factors.

## Usage

``` r
fit_glm(data, formula, peptide_id = NULL)
```

## Arguments

- data:

  Data frame containing the peptide data

- formula:

  A formula specifying the model

- peptide_id:

  Character, the peptide being analyzed (for diagnostics)

## Value

A list with components:

- converged:

  Logical, whether the model converged

- model:

  The fitted glm object (NULL if failed)

- coefficients:

  Model coefficients (NULL if failed)

- deviance:

  Residual deviance (NA if failed)

- peptide:

  The peptide ID
