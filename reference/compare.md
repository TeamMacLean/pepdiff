# Compare peptide abundances between conditions

Performs differential abundance analysis on proteomics data. Supports
three methods: GLM (default), ART (non-parametric), and pairwise tests.

## Usage

``` r
compare(data, ...)

# S3 method for class 'pepdiff_data'
compare(
  data,
  compare,
  ref,
  within = NULL,
  method = c("glm", "art", "pairwise"),
  test = c("wilcoxon", "bootstrap_t", "bayes_t", "rankprod"),
  alpha = 0.05,
  fdr_method = "BH",
  bf_threshold = 3,
  ...
)
```

## Arguments

- data:

  A pepdiff_data object from \[read_pepdiff()\]

- ...:

  Additional arguments passed to methods

- compare:

  Factor to compare (character string)

- ref:

  Reference level for comparisons

- within:

  Optional factor(s) to stratify by

- method:

  Analysis method: "glm" (default), "art", or "pairwise"

- test:

  For pairwise method: "wilcoxon", "bootstrap_t", "bayes_t", or
  "rankprod"

- alpha:

  Significance threshold (default 0.05). Used for p-value based tests.

- fdr_method:

  FDR correction method (default "BH"). Not applied for bayes_t.

- bf_threshold:

  Bayes factor threshold for significance (default 3). Only used when
  test = "bayes_t". BF \> threshold marks peptide as significant.

## Value

A pepdiff_results object containing:

- results:

  Tibble with peptide, gene_id, comparison, fold_change, log2_fc,
  p_value, fdr, significant. For bayes_t: p_value/fdr are NA, includes
  bf and evidence columns.

- comparisons:

  Tibble defining the comparisons made

- method:

  Statistical method used

- diagnostics:

  Model convergence information (for GLM/ART)

- params:

  Analysis parameters

- data:

  The original pepdiff_data object

- call:

  The function call

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple comparison
results <- compare(data, compare = "treatment", ref = "ctrl")

# Stratified comparison
results <- compare(data, compare = "treatment", ref = "ctrl", within = "timepoint")

# Pairwise test
results <- compare(data, compare = "treatment", ref = "ctrl",
                   method = "pairwise", test = "wilcoxon")

# Bayes factor test (uses bf_threshold instead of alpha/FDR)
results <- compare(data, compare = "treatment", ref = "ctrl",
                   method = "pairwise", test = "bayes_t", bf_threshold = 10)
} # }
```
