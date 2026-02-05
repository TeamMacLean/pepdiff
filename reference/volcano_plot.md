# volcano plot the data

draws a plot of peptide count against log fc at either protein or
peptide level for samples

## Usage

``` r
volcano_plot(
  l,
  log = FALSE,
  base = 2,
  sig_level = 0.05,
  metric = "bootstrap_t_p_val",
  option = "E",
  direction = -1
)
```

## Arguments

- l:

  list of results data frames, typically from \`compare_many()\`

- log:

  log the data

- base:

  base for logging

- sig_level:

  significance cutoff for colour

- metric:

  metric to use for significance

- option:

  viridis colour scheme key to use

- direction:

  viridis colour scheme direction (1/-1)

## Value

ggplot2 plot
