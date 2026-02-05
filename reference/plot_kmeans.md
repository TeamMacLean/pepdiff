# K-means cluster the data on the samples

Performs and draws a K-means cluster on the samples. Estimates number of
clusters as the product of the number of treatments and seconds. So
tries to group the bio reps together

## Usage

``` r
plot_kmeans(df, nstart = 25, iter.max = 1000)
```

## Arguments

- df:

  dataframe, typically from \`import_data()\`

- nstart:

  nstart points for \`kmeans()\` function

- iter.max:

  max iterations to perform for \`kmeans()\` function

## Value

ggplot2 plot
