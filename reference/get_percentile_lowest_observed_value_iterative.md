# get p values for contrast using normal percentile

get p values for contrast using normal percentile

Get p-values using normal percentile (legacy)

## Usage

``` r
get_percentile_lowest_observed_value_iterative(
  treatment,
  control,
  iters = 1000
)

get_percentile_lowest_observed_value_iterative(
  treatment,
  control,
  iters = 1000
)
```

## Arguments

- treatment:

  Treatment matrix

- control:

  Control matrix

- iters:

  Number of iterations

## Value

dataframe with one column \`norm_quantile_pval\`

Data frame with p-values
