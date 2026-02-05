# Extract contrasts from a fitted ART model

Uses art.con() from ARTool to compute contrasts.

## Usage

``` r
extract_contrasts_art(model, specs, adjust = "none")
```

## Arguments

- model:

  A fitted art object

- specs:

  Factor name or formula for contrasts

- adjust:

  P-value adjustment method

## Value

A tibble with contrast results
