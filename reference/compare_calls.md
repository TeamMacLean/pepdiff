# compare sets of significant peptides called by the used data

produces an UpSet plot showing intersections and set-size of the
different sets of significant peptides called by the methods used in the
provided result dataframe

## Usage

``` r
compare_calls(r, sig = 0.05)
```

## Arguments

- r:

  result dataframe typically from \`compare()\`

- sig:

  significance cut-off to select peptides

## Value

UpSet plot
