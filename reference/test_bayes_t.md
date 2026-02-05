# Bayes factor t-test for two groups

Computes a Bayes factor comparing the alternative hypothesis (group
difference) to the null hypothesis (no difference) using the JZS
(Jeffreys-Zellner-Siow) prior. Uses an analytical approximation for
computational efficiency.

## Usage

``` r
test_bayes_t(control, treatment, r_scale = 0.707)
```

## Arguments

- control:

  Numeric vector of control group values

- treatment:

  Numeric vector of treatment group values

- r_scale:

  Scale parameter for the Cauchy prior on effect size (default 0.707)

## Value

A list with components:

- bf:

  Bayes factor (BF10) - evidence for alternative vs null

- effect_size:

  Cohen's d effect size

- method:

  "bayes_t"

## Details

The Bayes factor is interpreted as: - BF10 \> 10: Strong evidence for
difference - BF10 \> 3: Moderate evidence for difference - BF10 0.33-3:
Inconclusive - BF10 \< 0.33: Moderate evidence for no difference - BF10
\< 0.1: Strong evidence for no difference

Unlike p-values, Bayes factors are NOT converted to pseudo-p-values. Use
\[classify_bf_evidence()\] to interpret BF values categorically.

## Examples

``` r
ctrl <- c(100, 120, 110, 105)
trt <- c(200, 220, 180, 210)
test_bayes_t(ctrl, trt)
#> $bf
#> [1] 1e+10
#> 
#> $effect_size
#> [1] 6.943651
#> 
#> $method
#> [1] "bayes_t"
#> 
```
