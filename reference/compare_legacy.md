# Legacy compare function (deprecated)

For each peptide this function carries out the selected tests to
determine peptides that are differentially abundant in the two
experiments specified. Before tests are performed the technical
replicates are merged and lowest observed value replacement for each
missing biological replicate is done. All comparisins are performed as
\`treatment / control\`

## Usage

``` r
compare_legacy(
  df,
  iters = 1000,
  treatment = NA,
  t_seconds = NA,
  control = NA,
  c_seconds = NA,
  tests = c("bootstrap_t")
)
```

## Arguments

- df:

  dataframe. Typically from \`import_data()\`

- iters:

  number of iterations to perform for iterative tests

- treatment:

  name of the experimental treatment to use as the \`treatment\`
  condition

- t_seconds:

  time point of the \`treatment\` condition to use

- control:

  name of the experimental treatment to use as the \`control\` condition

- c_seconds:

  time point of the \`control\` condition to use

- tests:

  character vector of tests to use, one or more of: \`norm_quantile\`,
  \`bootstrap_t\`, \`wilcoxon\`, \`kruskal-wallis\`, \`rank_product\`

## Value

dataframe with original and replaced quantification values, natural fold
change, biological replicates and p-value / fdr for each
