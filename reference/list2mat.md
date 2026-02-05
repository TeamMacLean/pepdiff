# converts a results object to a matrix as if for direct use in external heatmap functions

converts a results object to a matrix as if for direct use in external
heatmap functions

## Usage

``` r
list2mat(r, column = "fold_change")
```

## Arguments

- r:

  results object, usually from \`compare_many()\`

- column:

  column from results data to put into matrix, default = "fold_change"
