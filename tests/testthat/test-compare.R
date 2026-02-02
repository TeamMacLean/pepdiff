# Tests for compare.R - main analysis function

# =============================================================================
# Test Input Validation
# =============================================================================

test_that("compare requires pepdiff_data object", {
  expect_error(
    compare(data.frame(a = 1), compare = "treatment", ref = "ctrl"),
    class = "simpleError"
  )
})


test_that("compare requires valid compare factor", {
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

  expect_error(
    compare(imported, compare = "nonexistent", ref = "ctrl"),
    "not found"
  )
})


test_that("compare requires valid reference level", {
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

  expect_error(
    compare(imported, compare = "treatment", ref = "nonexistent"),
    "not found"
  )
})


test_that("compare uses first level as reference if not specified", {
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

  # Should work with message about default reference
  expect_message(
    compare(imported, compare = "treatment", method = "pairwise"),
    "reference level"
  )
})


# =============================================================================
# Test Pairwise Method
# =============================================================================

test_that("compare with pairwise wilcoxon works", {
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

  expect_s3_class(result, "pepdiff_results")
  expect_equal(result$method, "pairwise")
  expect_true("p_value" %in% names(result$results))
  expect_true("fdr" %in% names(result$results))
  expect_true("significant" %in% names(result$results))
})


test_that("compare with pairwise bootstrap_t works", {
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

  expect_s3_class(result, "pepdiff_results")
  expect_equal(result$params$test, "bootstrap_t")
})


test_that("compare with pairwise bayes_t works", {
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

  expect_s3_class(result, "pepdiff_results")
  expect_equal(result$params$test, "bayes_t")
})


test_that("compare with pairwise rankprod works", {
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
    test = "rankprod"
  )

  expect_s3_class(result, "pepdiff_results")
  expect_equal(result$params$test, "rankprod")
})


# =============================================================================
# Test GLM Method
# =============================================================================

test_that("compare with glm method works", {
  skip_if_not_installed("emmeans")

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
    method = "glm"
  )

  expect_s3_class(result, "pepdiff_results")
  expect_equal(result$method, "glm")
  expect_true("diagnostics" %in% names(result))
})


test_that("compare with glm includes convergence info", {
  skip_if_not_installed("emmeans")

  test_data <- make_test_data(n_peptides = 10, n_reps = 3, seed = 42)
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
    method = "glm"
  )

  expect_true("converged" %in% names(result$diagnostics))
})


# =============================================================================
# Test ART Method
# =============================================================================

test_that("compare with art method works", {
  skip_if_not_installed("ARTool")

  test_data <- make_test_data(n_peptides = 5, n_reps = 4, seed = 42)
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
    method = "art"
  )

  expect_s3_class(result, "pepdiff_results")
  expect_equal(result$method, "art")
})


# =============================================================================
# Test Results Structure
# =============================================================================

test_that("compare results have correct structure", {
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

  # Check results tibble structure
  expect_true("peptide" %in% names(result$results))
  expect_true("gene_id" %in% names(result$results))
  expect_true("comparison" %in% names(result$results))
  expect_true("fold_change" %in% names(result$results))
  expect_true("log2_fc" %in% names(result$results))
  expect_true("p_value" %in% names(result$results))
  expect_true("fdr" %in% names(result$results))
  expect_true("significant" %in% names(result$results))

  # Check comparisons tibble
  expect_true("comparison" %in% names(result$comparisons))

  # Check params
  expect_equal(result$params$compare, "treatment")
  expect_equal(result$params$ref, "ctrl")
  expect_equal(result$params$alpha, 0.05)
})


test_that("compare applies FDR correction within comparisons", {
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
    alpha = 0.05,
    fdr_method = "BH"
  )

  # FDR values should be >= p-values
  valid_idx <- !is.na(result$results$p_value) & !is.na(result$results$fdr)
  if (any(valid_idx)) {
    expect_true(all(result$results$fdr[valid_idx] >= result$results$p_value[valid_idx]))
  }
})


# =============================================================================
# Test with Factorial Design
# =============================================================================

test_that("compare works with factorial design", {
  skip_if_not_installed("emmeans")

  test_data <- make_factorial_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint"),
    replicate = "bio_rep"
  )

  # GLM with factorial design
  result <- compare(
    imported,
    compare = "treatment",
    ref = "ctrl",
    method = "glm"
  )

  expect_s3_class(result, "pepdiff_results")
})


# =============================================================================
# Test Output Methods
# =============================================================================

test_that("compare results can be printed", {
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
    method = "pairwise"
  )

  output <- capture.output(print(result))
  expect_true(length(output) > 0)
})


test_that("compare results can be summarized", {
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
    method = "pairwise"
  )

  output <- capture.output(summ <- summary(result))
  expect_true(length(output) > 0)
  expect_type(summ, "list")
})


# =============================================================================
# Test significant() accessor
# =============================================================================

test_that("significant() extracts significant results from compare output", {
  test_data <- make_test_data(
    n_peptides = 20,
    n_reps = 3,
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

  result <- compare(
    imported,
    compare = "treatment",
    ref = "ctrl",
    method = "pairwise",
    test = "wilcoxon"
  )

  sig <- significant(result)

  # All returned results should be significant
  expect_true(all(sig$fdr < 0.05 | is.na(sig$fdr)))
})
