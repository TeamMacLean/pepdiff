# get p values for contrast using Wilcoxon test

get p values for contrast using Wilcoxon test

Legacy: Get Wilcoxon test p-values for matrix data

## Usage

``` r
get_wilcoxon_percentile(treatment, control)

get_wilcoxon_percentile(treatment, control)
```

## Arguments

- treatment:

  Matrix of treatment data

- control:

  Matrix of control data

## Value

dataframe with two columns \`wilcoxon_p_val\` and \`wilcoxon_fdr\`

Data frame with wilcoxon_p_val and wilcoxon_fdr columns
