# works out if a peptide has at least one significant value across the experiment Composes a matrix of the \`metric\` significance values with peptides in rows, experiments in columns and works out if each peptide row has a value below the stated cut off

\#' returns a logical vector of length equal to row number of matrix

## Usage

``` r
get_sig_rows(l, metric = "bootstrap_t_pval", sig_level = 0.05)
```

## Arguments

- l:

  list of results, usually from \`compare_many()\`

- metric:

  the test metric used to determine significance one of:

- sig_level:

  significance level cutoff
