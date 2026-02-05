# Rank products test for two groups

\`r lifecycle::badge("deprecated")\`

This function is deprecated because Rank Products requires ranking
across ALL peptides, not within single peptides. The per-peptide
permutation approach produces unreliable p-values.

Use \`compare()\` with \`test = "rankprod"\` instead, which properly
uses the RankProd package to rank across all peptides.

## Usage

``` r
test_rankprod(control, treatment, n_perm = 1000, seed = NULL)
```

## Arguments

- control:

  Numeric vector of control group values

- treatment:

  Numeric vector of treatment group values

- n_perm:

  Number of permutations for p-value estimation (default 1000)

- seed:

  Optional random seed for reproducibility

## Value

A list with components:

- p_value_up:

  P-value for upregulation (treatment \> control)

- p_value_down:

  P-value for downregulation (treatment \< control)

- p_value:

  Combined two-sided p-value (minimum of up/down)

- rp_up:

  Rank product for upregulation

- rp_down:

  Rank product for downregulation

- method:

  "rankprod"

## Examples

``` r
ctrl <- c(100, 120, 110, 105)
trt <- c(200, 220, 180, 210)
if (FALSE) { # \dontrun{
test_rankprod(ctrl, trt, n_perm = 100)  # Deprecated
} # }
```
