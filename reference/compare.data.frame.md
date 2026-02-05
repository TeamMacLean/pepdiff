# Default compare method for legacy data frames

This method handles calls to compare() with data frames from the old
import_data() function. It issues a deprecation warning and delegates to
compare_legacy().

## Usage

``` r
# S3 method for class 'data.frame'
compare(data, ...)
```

## Arguments

- data:

  A data frame from import_data()

- ...:

  Arguments passed to compare_legacy()

## Value

Results from compare_legacy()
