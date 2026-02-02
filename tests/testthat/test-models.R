# Tests for models.R - GLM and ART model fitting

# =============================================================================
# Test fit_glm
# =============================================================================

test_that("fit_glm fits a simple model successfully", {
  data <- data.frame(
    value = c(100, 110, 105, 200, 210, 205),
    treatment = factor(c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
  )

  result <- fit_glm(data, value ~ treatment, peptide_id = "PEP_001")

  expect_true(result$converged)
  expect_false(is.null(result$model))
  expect_equal(result$peptide, "PEP_001")
})


test_that("fit_glm fails with zero values", {
  data <- data.frame(
    value = c(0, 110, 105, 200, 210, 205),
    treatment = factor(c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
  )

  result <- fit_glm(data, value ~ treatment, peptide_id = "PEP_001")

  expect_false(result$converged)
  expect_true(grepl("positive", result$error, ignore.case = TRUE))
})


test_that("fit_glm fails with negative values", {
  data <- data.frame(
    value = c(-10, 110, 105, 200, 210, 205),
    treatment = factor(c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
  )

  result <- fit_glm(data, value ~ treatment, peptide_id = "PEP_001")

  expect_false(result$converged)
})


test_that("fit_glm handles factorial design", {
  data <- data.frame(
    value = c(100, 110, 200, 210, 150, 160, 300, 310),
    treatment = factor(rep(c("ctrl", "trt"), each = 4)),
    timepoint = factor(rep(c("0h", "24h"), 4))
  )

  result <- fit_glm(data, value ~ treatment * timepoint, peptide_id = "PEP_001")

  expect_true(result$converged)
  expect_true(length(result$coefficients) > 1)
})


# =============================================================================
# Test extract_contrasts_glm
# =============================================================================

test_that("extract_contrasts_glm extracts pairwise contrasts", {
  skip_if_not_installed("emmeans")

  data <- data.frame(
    value = c(100, 110, 105, 200, 210, 205),
    treatment = factor(c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
  )

  model <- stats::glm(value ~ treatment, data = data, family = stats::Gamma(link = "log"))
  result <- extract_contrasts_glm(model, specs = "treatment")

  expect_s3_class(result, "tbl_df")
  expect_true("contrast" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("fold_change" %in% names(result))
})


test_that("extract_contrasts_glm handles reference level", {
  skip_if_not_installed("emmeans")

  data <- data.frame(
    value = c(100, 110, 105, 200, 210, 205, 150, 160, 155),
    treatment = factor(c(rep("ctrl", 3), rep("trt1", 3), rep("trt2", 3)))
  )

  model <- stats::glm(value ~ treatment, data = data, family = stats::Gamma(link = "log"))

  # Treatment vs control contrasts
  result <- extract_contrasts_glm(model, specs = "treatment", ref = 1)

  expect_equal(nrow(result), 2)  # trt1 vs ctrl, trt2 vs ctrl
})


# =============================================================================
# Test fit_and_extract_glm
# =============================================================================

test_that("fit_and_extract_glm works end-to-end", {
  skip_if_not_installed("emmeans")

  data <- data.frame(
    value = c(100, 110, 105, 200, 210, 205),
    treatment = factor(c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
  )

  result <- fit_and_extract_glm(
    data = data,
    response = "value",
    factors = "treatment",
    compare = "treatment",
    ref = 1,
    peptide_id = "PEP_001"
  )

  expect_equal(result$peptide, "PEP_001")
  expect_true(result$converged)
  expect_s3_class(result$contrasts, "tbl_df")
})


test_that("fit_and_extract_glm handles failure gracefully", {
  data <- data.frame(
    value = c(0, 0, 0, 200, 210, 205),
    treatment = factor(c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
  )

  result <- fit_and_extract_glm(
    data = data,
    response = "value",
    factors = "treatment",
    compare = "treatment",
    ref = 1,
    peptide_id = "PEP_001"
  )

  expect_false(result$converged)
  expect_null(result$contrasts)
})


# =============================================================================
# Test fit_art (if ARTool available)
# =============================================================================

test_that("fit_art works when ARTool is installed", {
  skip_if_not_installed("ARTool")

  data <- data.frame(
    value = c(100, 110, 105, 108, 200, 210, 205, 208),
    treatment = factor(rep(c("ctrl", "trt"), each = 4)),
    subject = factor(rep(1:4, 2))
  )

  # ART needs Error() term for repeated measures, use simpler formula
  result <- fit_art(data, value ~ treatment, peptide_id = "PEP_001")

  expect_true(result$converged)
  expect_false(is.null(result$model))
})


test_that("fit_art fails gracefully without ARTool", {
  skip_if(requireNamespace("ARTool", quietly = TRUE), "ARTool is installed")

  data <- data.frame(
    value = c(100, 110, 200, 210),
    treatment = factor(c("ctrl", "ctrl", "trt", "trt"))
  )

  result <- fit_art(data, value ~ treatment, peptide_id = "PEP_001")

  expect_false(result$converged)
  expect_true(grepl("ARTool", result$error))
})


# =============================================================================
# Test run_models
# =============================================================================

test_that("run_models processes all peptides", {
  skip_if_not_installed("emmeans")

  # Create test data
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

  results <- run_models(imported, compare = "treatment", ref = 1, method = "glm")

  expect_equal(length(results), length(imported$peptides))
  expect_true(all(names(results) %in% imported$peptides))
})


test_that("run_models handles mixed convergence", {
  skip_if_not_installed("emmeans")

  # Create data where some peptides might have issues
  test_data <- make_test_data(n_peptides = 5, n_reps = 3, na_rate = 0.2, seed = 42)
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

  results <- run_models(imported, compare = "treatment", ref = 1, method = "glm")

  # Check that all peptides have results (even if some didn't converge)
  expect_equal(length(results), length(imported$peptides))

  # At least some should have converged
  n_converged <- sum(vapply(results, function(x) x$converged, logical(1)))
  expect_true(n_converged > 0)
})


# =============================================================================
# Test extract_diagnostics
# =============================================================================

test_that("extract_diagnostics creates summary tibble", {
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

  model_results <- run_models(imported, compare = "treatment", ref = 1, method = "glm")
  diagnostics <- extract_diagnostics(model_results)

  expect_s3_class(diagnostics, "tbl_df")
  expect_true("peptide" %in% names(diagnostics))
  expect_true("converged" %in% names(diagnostics))
  expect_equal(nrow(diagnostics), length(imported$peptides))
})


# =============================================================================
# Test helper functions
# =============================================================================

test_that("has_sufficient_variation detects constant values", {
  expect_false(has_sufficient_variation(c(100, 100, 100)))
  expect_true(has_sufficient_variation(c(100, 200, 150)))
})


test_that("has_sufficient_variation handles NA", {
  expect_true(has_sufficient_variation(c(100, NA, 200, 150)))
  expect_false(has_sufficient_variation(c(NA, NA)))
})


test_that("has_sufficient_variation needs minimum data", {
  expect_false(has_sufficient_variation(c(100)))
  expect_false(has_sufficient_variation(c(100, 200)))
  expect_true(has_sufficient_variation(c(100, 200, 150)))
})
