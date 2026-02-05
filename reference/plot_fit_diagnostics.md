# Plot GLM fit diagnostics

Creates a multi-panel diagnostic plot to help assess whether GLM models
fit the data well. This is useful for deciding whether to use GLM or
switch to ART (Aligned Rank Transform).

## Usage

``` r
plot_fit_diagnostics(
  results,
  n_sample = 6,
  deviance_threshold = NULL,
  full_qq = FALSE
)
```

## Arguments

- results:

  A \`pepdiff_results\` object from \`compare()\` with \`method =
  "glm"\`

- n_sample:

  Number of peptides to show in sample residual plots (default 6)

- deviance_threshold:

  Optional threshold for flagging high-deviance peptides. If NULL
  (default), uses the 95th percentile of deviance values.

- full_qq:

  Deprecated. Residuals are now stored during \`compare()\` so accurate
  QQ plots are always available without refitting.

## Value

Invisibly returns a list with:

- plot:

  The diagnostic plot (ggplot/cowplot grid)

- flagged:

  Tibble of peptides with potential fit issues

- summary:

  List with summary statistics (n_flagged, median_deviance, etc.)

## Details

The function generates a 4-panel diagnostic plot:

\*\*Panel 1: Deviance Distribution\*\* - Histogram showing the
distribution of residual deviance across all converged peptides. A long
right tail suggests some peptides fit poorly.

\*\*Panel 2: Deviance vs Fold Change\*\* - Scatter plot of deviance
against absolute log2 fold change. If high-deviance points cluster at
extreme fold changes, this may indicate outlier-driven "significant"
results.

\*\*Panel 3: Sample Residual Plots\*\* - Residuals vs fitted values for
a sample of peptides (2 with highest deviance, 2 median, 2 lowest). Look
for random scatter around zero; patterns or funnels indicate poor fit.

\*\*Panel 4: Pooled QQ Plot\*\* - Quantile-quantile plot of pooled
residuals. Points should fall on the diagonal line. S-curves indicate
heavy tails (consider ART), systematic deviation suggests wrong
distributional assumption.

## Interpretation

\*\*Use GLM when:\*\* - Deviance distribution looks reasonable (few
flagged peptides) - No systematic patterns in residual plots - QQ plot
is reasonably linear

\*\*Consider ART when:\*\* - Many peptides (\>15 - Residual plots show
systematic curves or funnels - QQ plot shows heavy tails (S-curve)

## See also

\[compare()\] for running the analysis, \`vignette("checking_fit")\` for
detailed guidance on interpreting diagnostics

## Examples

``` r
if (FALSE) { # \dontrun{
# Run GLM analysis
results <- compare(dat, compare = "treatment", ref = "ctrl", method = "glm")

# Check fit diagnostics
diag <- plot_fit_diagnostics(results)

# View flagged peptides
diag$flagged

# Get summary statistics
diag$summary
} # }
```
