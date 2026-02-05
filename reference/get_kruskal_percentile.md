# get p values for contrast using Kruskal-Wallis test

get p values for contrast using Kruskal-Wallis test

Legacy: Get Kruskal-Wallis test p-values for matrix data

## Usage

``` r
get_kruskal_percentile(treatment, control)

get_kruskal_percentile(treatment, control)
```

## Arguments

- treatment:

  Matrix of treatment data

- control:

  Matrix of control data

## Value

dataframe with two columns \`kruskal_p_val\` and \`kruskal_fdr\`

Data frame with kruskal_p_val and kruskal_fdr columns
