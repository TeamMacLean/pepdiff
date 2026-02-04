# GLM/ART Stratification Fix Specification

## Problem

The `within` parameter is accepted by `compare()` but silently ignored for GLM and ART methods. The helper function `compare_within_strata()` exists and is complete, but is never called from the main dispatch.

## Current Behavior (broken)

```r
compare(dat, compare = "treatment", ref = "ctrl", within = "timepoint", method = "glm")
# Returns: main effect only (averaged across timepoints)
# The "within" parameter is ignored
```

## Expected Behavior

```r
compare(dat, compare = "treatment", ref = "ctrl", within = "timepoint", method = "glm")
# Returns: treatment effect at 0h AND treatment effect at 24h (separate rows per peptide)
```

## Results Structure

When `within = "timepoint"`, each peptide gets one row per stratum:

| peptide | comparison | timepoint | fold_change | p_value | fdr | significant |
|---------|------------|-----------|-------------|---------|-----|-------------|
| PEP_001 | trt vs ctrl | 0h | 2.9 | 0.001 | 0.01 | TRUE |
| PEP_001 | trt vs ctrl | 24h | 3.1 | 0.002 | 0.02 | TRUE |
| PEP_002 | trt vs ctrl | 0h | 1.1 | 0.45 | 0.50 | FALSE |
| PEP_002 | trt vs ctrl | 24h | 0.9 | 0.55 | 0.60 | FALSE |

## Design Decisions

1. **FDR correction**: Applied within each stratum separately (not globally)
2. **Comparison column**: Stays as "trt vs ctrl" - stratum info in separate column(s)
3. **Multiple within factors**: Supported - `within = c("timepoint", "genotype")` creates one row per combination

## The Fix

In `compare.pepdiff_data()` (R/compare.R ~line 106), change the dispatch logic:

```r
# Current (broken):
if (method == "pairwise") {
  results <- compare_pairwise(...)
} else if (method == "glm") {
  results <- compare_glm(...)  # ignores within
} else if (method == "art") {
  results <- compare_art(...)  # ignores within
}

# Fixed:
if (!is.null(within)) {
  # Stratified analysis - works for all methods
  results <- compare_within_strata(data, compare, ref, within, method, test, alpha, fdr_method)
} else if (method == "pairwise") {
  results <- compare_pairwise(...)
} else if (method == "glm") {
  results <- compare_glm(...)
} else if (method == "art") {
  results <- compare_art(...)
}
```

## Additional Changes Needed

1. **compare_within_strata()**: Add `bf_threshold` parameter for bayes_t pairwise tests
2. **comparisons tibble**: Should reflect strata when `within` is used
3. **print method**: Should show stratum information

## Test Cases

### Test 1: GLM stratification produces separate results per stratum

```r
test_that("GLM with within produces results for each stratum", {
  # Setup: 2x2 factorial data
  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "glm")

  # Should have 2 rows per peptide (one per timepoint)
  expect_equal(nrow(results$results), 2 * length(dat$peptides))

  # Should have timepoint column

  expect_true("timepoint" %in% names(results$results))

  # Should have both timepoint levels
  expect_setequal(unique(results$results$timepoint), c("0h", "24h"))
})
```

### Test 2: Stratified results detect interaction correctly

```r
test_that("GLM stratification detects interaction pattern", {
  # Using our vignette data with Group C (interaction: effect at 24h only)
  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "glm")

  group_c <- paste0("PEP_", sprintf("%03d", 21:30))

  group_c_results <- results$results %>%
    filter(peptide %in% group_c) %>%
    group_by(timepoint) %>%
    summarise(median_fc = median(fold_change))

  # At 0h: FC should be ~1 (no effect)
  fc_0h <- group_c_results$median_fc[group_c_results$timepoint == "0h"]
  expect_true(fc_0h > 0.7 && fc_0h < 1.5)

  # At 24h: FC should be ~4 (the true effect)
  fc_24h <- group_c_results$median_fc[group_c_results$timepoint == "24h"]
  expect_true(fc_24h > 3 && fc_24h < 5)
})
```

### Test 3: FDR applied within stratum

```r
test_that("FDR correction is applied within each stratum", {
  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "glm")

  # Check FDR is computed separately per stratum
  for (tp in c("0h", "24h")) {
    stratum_results <- results$results %>% filter(timepoint == tp)
    expected_fdr <- p.adjust(stratum_results$p_value, method = "BH")
    expect_equal(stratum_results$fdr, expected_fdr)
  }
})
```

### Test 4: Multiple within factors

```r
test_that("Multiple within factors work", {
  # Would need 3-factor data for this test
  # within = c("timepoint", "genotype") should produce
  # one row per peptide per timepoint-genotype combination
})
```

### Test 5: Pairwise with within still works

```r
test_that("Pairwise method with within still works", {
  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "pairwise", test = "wilcoxon")

  expect_true("timepoint" %in% names(results$results))
  expect_equal(nrow(results$results), 2 * length(dat$peptides))
})
```

## Files to Modify

1. `R/compare.R` - Add dispatch to `compare_within_strata()` when `within` is not NULL
2. `R/compare.R` - Update `compare_within_strata()` to handle `bf_threshold`
3. `R/results.R` - Update print method to show stratum info (if needed)
4. `tests/testthat/test-stratification.R` - New test file with above tests

## Success Criteria

1. `compare(..., within = "timepoint", method = "glm")` returns stratified results
2. Group C peptides show FC ≈ 1 at 0h, FC ≈ 4 at 24h (not diluted to ~2)
3. FDR applied within each stratum
4. All existing tests still pass
5. `devtools::check()` has 0 errors/warnings
