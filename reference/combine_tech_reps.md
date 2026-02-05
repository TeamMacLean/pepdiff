# Combine technical replicates

Explicitly combines technical replicates by averaging values within each
combination of peptide, factors, and biological replicate.

## Usage

``` r
combine_tech_reps(data, fun = mean)
```

## Arguments

- data:

  A pepdiff_data object with a tech_rep column

- fun:

  Function to use for combining (default: mean)

## Value

A pepdiff_data object with technical replicates combined

## Examples

``` r
if (FALSE) { # \dontrun{
# Import data with technical replicates
dat <- read_pepdiff(..., tech_rep = "tech_rep")

# Combine by averaging
dat <- combine_tech_reps(dat)

# Or combine by taking median
dat <- combine_tech_reps(dat, fun = median)
} # }
```
