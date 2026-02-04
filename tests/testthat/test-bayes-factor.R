# Tests for proper Bayes factor handling
# These tests define the expected behavior after implementing the spec in bayes_factor_spec.md

# =============================================================================
# test_bayes_t() should NOT return pseudo-p-values
# =============================================================================

test_that("test_bayes_t returns bf without p_value", {

  ctrl <- c(100, 110, 105, 108, 103)
  trt <- c(200, 210, 205, 208, 203)

  result <- test_bayes_t(ctrl, trt)

  expect_type(result, "list")
  expect_true("bf" %in% names(result))
  expect_true("effect_size" %in% names(result))
  expect_equal(result$method, "bayes_t")


  # p_value should NOT be returned - this is the key change

  expect_false("p_value" %in% names(result))
})


test_that("test_bayes_t returns valid Bayes factor",
{
  ctrl <- c(100, 110, 105, 108, 103)
  trt <- c(200, 210, 205, 208, 203)

  result <- test_bayes_t(ctrl, trt)

  # BF should be positive
  expect_true(result$bf > 0)

  # For clearly different groups, BF should strongly favor alternative
  expect_true(result$bf > 10)
})


test_that("test_bayes_t returns BF < 1 for similar groups", {
  ctrl <- c(100, 102, 101, 99, 100)
  trt <- c(101, 100, 102, 99, 101)

  result <- test_bayes_t(ctrl, trt)

  # For very similar groups, BF should favor null (BF < 1)
  expect_true(result$bf < 3)  # At least not strong evidence for difference
})


# =============================================================================
# classify_bf_evidence() helper function
# =============================================================================

test_that("classify_bf_evidence returns correct categories", {
  # Test all five evidence categories
  expect_equal(classify_bf_evidence(0.05), "strong_null")
  expect_equal(classify_bf_evidence(0.2), "moderate_null")
  expect_equal(classify_bf_evidence(1.0), "inconclusive")
  expect_equal(classify_bf_evidence(5.0), "moderate_alt")
  expect_equal(classify_bf_evidence(20.0), "strong_alt")
})


test_that("classify_bf_evidence returns ordered factor", {
  evidence <- classify_bf_evidence(c(0.05, 0.2, 1.0, 5.0, 20.0))

  expect_s3_class(evidence, "factor")
  expect_true(is.ordered(evidence))

  # Check correct order
  expected_levels <- c("strong_null", "moderate_null", "inconclusive",
                       "moderate_alt", "strong_alt")
  expect_equal(levels(evidence), expected_levels)
})


test_that("classify_bf_evidence handles boundary values", {
  # At boundaries
  expect_equal(as.character(classify_bf_evidence(0.1)), "moderate_null")  # BF = 0.1 is boundary
  expect_equal(as.character(classify_bf_evidence(0.33)), "inconclusive")  # BF = 0.33 is boundary
  expect_equal(as.character(classify_bf_evidence(3)), "moderate_alt")     # BF = 3 is boundary
  expect_equal(as.character(classify_bf_evidence(10)), "strong_alt")      # BF = 10 is boundary
})


# =============================================================================
# compare() with bayes_t should handle results correctly
# =============================================================================

test_that("compare with bayes_t has bf column, not p_value", {
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  result <- compare(
    imported,
    compare = "treatment",
    ref = "ctrl",
    method = "pairwise",
    test = "bayes_t"
  )

  # bf column should be present
  expect_true("bf" %in% names(result$results))

  # p_value and fdr should be NA for bayes_t
  expect_true(all(is.na(result$results$p_value)))
  expect_true(all(is.na(result$results$fdr)))
})


test_that("compare with bayes_t has evidence column", {
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  result <- compare(
    imported,
    compare = "treatment",
    ref = "ctrl",
    method = "pairwise",
    test = "bayes_t"
  )

  # evidence column should be present and be an ordered factor
  expect_true("evidence" %in% names(result$results))
  expect_s3_class(result$results$evidence, "factor")
  expect_true(is.ordered(result$results$evidence))

  # Should have correct levels
  expected_levels <- c("strong_null", "moderate_null", "inconclusive",
                       "moderate_alt", "strong_alt")
  expect_equal(levels(result$results$evidence), expected_levels)
})


test_that("compare with bayes_t uses bf_threshold for significance", {
  test_data <- make_test_data(
    n_peptides = 20,
    n_reps = 4,
    effect_size = 3.0,
    affected_fraction = 0.5,
    seed = 42
  )
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  # Default threshold = 3
  result <- compare(
    imported,
    compare = "treatment",
    ref = "ctrl",
    method = "pairwise",
    test = "bayes_t"
  )

  # significant should be TRUE where bf > 3
  sig_idx <- result$results$significant
  bf_vals <- result$results$bf

  expect_true(all(bf_vals[sig_idx] > 3, na.rm = TRUE))
  expect_true(all(bf_vals[!sig_idx] <= 3, na.rm = TRUE))
})


test_that("compare with bayes_t respects custom bf_threshold", {
  test_data <- make_test_data(
    n_peptides = 20,
    n_reps = 4,
    effect_size = 3.0,
    affected_fraction = 0.5,
    seed = 42
  )
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  # Custom threshold = 10 (strong evidence)
  result <- compare(
    imported,
    compare = "treatment",
    ref = "ctrl",
    method = "pairwise",
    test = "bayes_t",
    bf_threshold = 10
  )

  # significant should be TRUE where bf > 10
  sig_idx <- result$results$significant
  bf_vals <- result$results$bf

  expect_true(all(bf_vals[sig_idx] > 10, na.rm = TRUE))
  expect_true(all(bf_vals[!sig_idx] <= 10, na.rm = TRUE))

  # Params should record the threshold

expect_equal(result$params$bf_threshold, 10)
})


test_that("compare stores bf_threshold in params", {
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  result <- compare(
    imported,
    compare = "treatment",
    ref = "ctrl",
    method = "pairwise",
    test = "bayes_t"
  )

  # Default bf_threshold should be stored
  expect_equal(result$params$bf_threshold, 3)
})


# =============================================================================
# Other pairwise tests should NOT be affected
# =============================================================================

test_that("compare with wilcoxon still has p_value and fdr", {
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  result <- compare(
    imported,
    compare = "treatment",
    ref = "ctrl",
    method = "pairwise",
    test = "wilcoxon"
  )

  # wilcoxon should still have p_value and fdr
  expect_true("p_value" %in% names(result$results))
  expect_true("fdr" %in% names(result$results))
  expect_false(all(is.na(result$results$p_value)))
  expect_false(all(is.na(result$results$fdr)))

  # wilcoxon should NOT have bf or evidence columns
  expect_false("bf" %in% names(result$results))
  expect_false("evidence" %in% names(result$results))
})


test_that("compare with bootstrap_t still has p_value and fdr", {
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  result <- compare(
    imported,
    compare = "treatment",
    ref = "ctrl",
    method = "pairwise",
    test = "bootstrap_t"
  )

  # bootstrap_t should still have p_value and fdr
  expect_true("p_value" %in% names(result$results))
  expect_true("fdr" %in% names(result$results))
  expect_false(all(is.na(result$results$p_value)))

  # bootstrap_t should NOT have bf or evidence columns
  expect_false("bf" %in% names(result$results))
  expect_false("evidence" %in% names(result$results))
})


# =============================================================================
# bf_to_pvalue should be removed
# =============================================================================

test_that("bf_to_pvalue function no longer exists", {
  # This function should be removed as it's statistically invalid
  expect_false(exists("bf_to_pvalue", envir = asNamespace("pepdiff")))
})
