# Perform kmeans of a dataset using just data in selected columns, then return matrices of all columns

Perform kmeans of a dataset using just data in selected columns, then
return matrices of all columns

## Usage

``` r
kmeans_by_selected_cols(
  l,
  cols = NULL,
  log = TRUE,
  base = 2,
  sig_only = TRUE,
  sig_level = 0.05,
  metric = "bootstrap_t_p_val",
  k = NA,
  nstart = 25,
  itermax = 1000
)
```

## Arguments

- l:

  list of results, usually from \`compare_many()\`

- cols:

  names of columns to perform the k-means with

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

- k:

  number of clusters to make

- nstart:

  nstart value for \`kmeans()\`

- itermax:

  number of \`kmeans()\` iterations (1000)

## Value

list of matrices
