# returns a matrix of fold change values

Computes the fold change relative to the control sample and returns a
matrix with comparisons in columns and peptides in rows. Use this if you
want data for a customised heatmap

## Usage

``` r
fold_change_matrix(
  l,
  log = TRUE,
  base = 2,
  sig_only = FALSE,
  sig_level = 0.05,
  metric = "bootstrap_t_fdr"
)
```

## Arguments

- l:

  list of results, usually from \`compare_many()\`

- log:

  whether to log the data

- base:

  base used in logging (default = 2)

- sig_only:

  return only rows with 1 or more values significant at \`sig_level\` of
  \`metric\`

- sig_level:

  significance level cutoff

- metric:

  the test metric used to determine significance one of:
  \`bootstrap_t_p_val\`, \`bootstrap_t_fdr\` \`wilcoxon_p_val\`,
  \`wilcoxon_fdr\` \`kruskal_p_val\`, \`kruskal_fdr\`
  \`rank_prod_p1_p_val\`, \`rank_prod_p2_p_val\`, \`rank_prod_p1_fdr\`,
  \`rank_prod_p2_fdr\`.

## Value

matrix
