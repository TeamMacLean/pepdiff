# Wilcoxon rank-sum test for two groups

Performs a two-sample Wilcoxon rank-sum test (Mann-Whitney U test) to
compare abundance values between control and treatment groups.

## Usage

``` r
test_wilcoxon(control, treatment, ...)
```

## Arguments

- control:

  Numeric vector of control group values

- treatment:

  Numeric vector of treatment group values

- ...:

  Additional arguments passed to \[stats::wilcox.test()\]

## Value

A list with components:

- p_value:

  The p-value from the test

- statistic:

  The test statistic W

- method:

  "wilcoxon"

## Examples

``` r
ctrl <- c(100, 120, 110, 105)
trt <- c(200, 220, 180, 210)
test_wilcoxon(ctrl, trt)
#> $p_value
#> [1] 0.02857143
#> 
#> $statistic
#> [1] 16
#> 
#> $method
#> [1] "wilcoxon"
#> 
```
