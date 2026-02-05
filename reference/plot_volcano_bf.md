# Volcano plot for Bayes factor results

Creates a volcano plot with log10(BF) on the y-axis instead of
-log10(p-value). Reference lines are drawn at BF thresholds (3, 10) and
their reciprocals (0.33, 0.1).

## Usage

``` r
plot_volcano_bf(results, comparison = NULL, bf_threshold = 3, fc_threshold = 1)
```

## Arguments

- results:

  A pepdiff_results object from bayes_t test

- comparison:

  Optional comparison name to filter by

- bf_threshold:

  BF threshold for coloring (default 3)

- fc_threshold:

  Fold change threshold for labeling (default 1)

## Value

A ggplot object
