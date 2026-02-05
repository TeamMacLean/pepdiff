# compare many combinations of treatment and control

for each combination of treatment and control condition, runs the
\`compare()\` function and collates the results

## Usage

``` r
compare_many(df, comparison, iters = 1000, tests = c("bootstrap_t"))
```

## Arguments

- df:

  dataframe. Typically from \`import_data()\`

- comparison:

  path to file or dataframe of comparisons with columns treatment,
  t_seconds, control, c_seconds

- iters:

  number of iterations to perform for iterative tests

- tests:

  character vector of tests to use, one or more of: \`norm_quantile\`,
  \`bootstrap_t\`, \`wilcoxon\`, \`kruskal-wallis\`, \`rank_product\`
