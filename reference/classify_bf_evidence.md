# Classify Bayes factor into evidence categories

Converts numeric Bayes factors into categorical evidence levels
following conventional thresholds (Jeffreys, 1961; Lee & Wagenmakers,
2013).

## Usage

``` r
classify_bf_evidence(bf)
```

## Arguments

- bf:

  Numeric vector of Bayes factors (BF10)

## Value

An ordered factor with levels:

- strong_null:

  BF \< 0.1 - Strong evidence for null hypothesis

- moderate_null:

  BF 0.1-0.33 - Moderate evidence for null

- inconclusive:

  BF 0.33-3 - Evidence is inconclusive

- moderate_alt:

  BF 3-10 - Moderate evidence for alternative

- strong_alt:

  BF \> 10 - Strong evidence for alternative

## Examples

``` r
classify_bf_evidence(c(0.05, 0.2, 1, 5, 20))
#> [1] strong_null   moderate_null inconclusive  moderate_alt  strong_alt   
#> 5 Levels: strong_null < moderate_null < inconclusive < ... < strong_alt
```
