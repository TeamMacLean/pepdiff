# draw density plots for data

Plot density of quantities in data for each treatment, seconds and
biological replicate

## Usage

``` r
plot_quant_distributions(df, log = FALSE, base = 2)
```

## Arguments

- df:

  dataframe; typically from \`import_data()\`

- log:

  perform log transform of data

- base:

  base to use in log transform

## Value

ggplot2 plot
