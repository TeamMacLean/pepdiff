# Fit an Aligned Rank Transform model for a single peptide

Fits a non-parametric ART model which is suitable for factorial designs
when parametric assumptions are not met.

## Usage

``` r
fit_art(data, formula, peptide_id = NULL)
```

## Arguments

- data:

  Data frame containing the peptide data

- formula:

  A formula specifying the model

- peptide_id:

  Character, the peptide being analyzed

## Value

A list with components:

- converged:

  Logical, whether fitting succeeded

- model:

  The fitted art object (NULL if failed)

- peptide:

  The peptide ID

- error:

  Error message if failed
