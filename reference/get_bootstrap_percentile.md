# get p values for contrast using boostrap t test

get p values for contrast using boostrap t test

Legacy: Get bootstrap t-test p-values for matrix data

## Usage

``` r
get_bootstrap_percentile(treatment, control, iters = 1000)

get_bootstrap_percentile(treatment, control, iters = 1000)
```

## Arguments

- treatment:

  Matrix of treatment data (rows = peptides, cols = replicates)

- control:

  Matrix of control data

- iters:

  Number of bootstrap iterations

## Value

dataframe with two columns \`bootstrap_t_p_val\` and \`bootstrap_t_fdr\`

Data frame with bootstrap_t_p_val and bootstrap_t_fdr columns
