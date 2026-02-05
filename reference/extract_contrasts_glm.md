# Extract contrasts from a fitted GLM using emmeans

Uses the emmeans package to compute estimated marginal means and
contrasts between factor levels.

## Usage

``` r
extract_contrasts_glm(model, specs, ref = NULL, adjust = "none")
```

## Arguments

- model:

  A fitted glm object

- specs:

  Character vector of factor names to compare, or emmeans formula

- ref:

  Reference level for pairwise comparisons

- adjust:

  P-value adjustment method (default "none" - FDR done later)

## Value

A tibble with columns:

- contrast:

  Description of the contrast

- estimate:

  Log fold change estimate

- se:

  Standard error

- z_ratio:

  Z statistic

- p_value:

  P-value

- fold_change:

  Back-transformed fold change
