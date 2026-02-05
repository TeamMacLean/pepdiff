# Calculate mean fold change for peptide

For each peptide calculate the mean quantification over all experiments
then get the natural scale fold change.

## Usage

``` r
mean_fold_change(treatment, control)
```

## Arguments

- treatment:

  Vector or matrix of treatment data

- control:

  Vector or matrix of control data

## Value

Vector of mean fold changes: mean(treatment) / mean(control)
