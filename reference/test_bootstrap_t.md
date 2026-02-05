# Bootstrap t-test for two groups

Performs a bootstrap-based t-test comparing two groups. This is more
robust than a standard t-test when assumptions of normality may not
hold.

## Usage

``` r
test_bootstrap_t(control, treatment, n_boot = 1000, seed = NULL)
```

## Arguments

- control:

  Numeric vector of control group values

- treatment:

  Numeric vector of treatment group values

- n_boot:

  Number of bootstrap iterations (default 1000)

- seed:

  Optional random seed for reproducibility

## Value

A list with components:

- p_value:

  The two-sided p-value

- t_obs:

  The observed t-statistic

- method:

  "bootstrap_t"

## Examples

``` r
ctrl <- c(100, 120, 110, 105)
trt <- c(200, 220, 180, 210)
test_bootstrap_t(ctrl, trt, n_boot = 500)
#> $p_value
#> [1] 0.006
#> 
#> $t_obs
#> [1] 9.819805
#> 
#> $method
#> [1] "bootstrap_t"
#> 
```
