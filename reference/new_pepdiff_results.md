# Create a new pepdiff_results object

Low-level constructor for pepdiff_results objects. Typically created by
\[compare()\] rather than directly.

## Usage

``` r
new_pepdiff_results(
  results,
  comparisons,
  method,
  diagnostics,
  params,
  data,
  call
)
```

## Arguments

- results:

  Tibble with results in long format

- comparisons:

  Tibble defining the comparisons made

- method:

  Character, the statistical method used

- diagnostics:

  Tibble with model diagnostics (convergence, etc.)

- params:

  List of parameters used in the analysis

- data:

  The original pepdiff_data object

- call:

  The original function call

## Value

A pepdiff_results object
