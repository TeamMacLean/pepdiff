# JZS Bayes factor approximation

Computes BF10 using the Savage-Dickey density ratio approximation.

## Usage

``` r
jzs_bf_approx(t_stat, n_eff, r_scale = 0.707)
```

## Arguments

- t_stat:

  T-statistic

- n_eff:

  Effective sample size

- r_scale:

  Scale parameter for Cauchy prior

## Value

Bayes factor (BF10)
