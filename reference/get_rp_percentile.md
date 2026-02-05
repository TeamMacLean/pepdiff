# get p values for contrast using Rank Products test

get p values for contrast using Rank Products test

Legacy: Get Rank Products test p-values for matrix data

## Usage

``` r
get_rp_percentile(treatment, control)

get_rp_percentile(treatment, control)
```

## Arguments

- treatment:

  Matrix of treatment data

- control:

  Matrix of control data

## Value

dataframe with four columns, two for the test each way from RankProducts
\`rank_prod_p1_p_val\`, \`rank_prod_p2_p_val\` and \`rank_prod_p1_fdr\`,
\`rank_prod_p2_fdr\`.

Data frame with rank product p-values and FDR
