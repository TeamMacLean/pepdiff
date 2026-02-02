# Tests for utils.R

test_that("validate_factors checks for missing factors", {
  data <- data.frame(
    treatment = c("ctrl", "ctrl", "trt", "trt"),
    value = c(1, 2, 3, 4)
  )

  # Should pass with existing factor
expect_true(validate_factors(data, "treatment"))

  # Should fail with non-existent factor
  expect_error(
    validate_factors(data, "nonexistent"),
    "Factor columns not found"
  )
  expect_error(
    validate_factors(data, c("treatment", "missing")),
    "missing"
  )
})

test_that("validate_factors requires at least 2 levels", {
  data <- data.frame(
    treatment = c("ctrl", "ctrl", "ctrl"),
    value = c(1, 2, 3)
  )

  expect_error(
    validate_factors(data, "treatment"),
    "at least 2 levels"
  )
})

test_that("validate_factors handles multiple factors", {
  data <- data.frame(
    treatment = c("ctrl", "ctrl", "trt", "trt"),
    timepoint = c("0h", "24h", "0h", "24h"),
    value = c(1, 2, 3, 4)
  )

  expect_true(validate_factors(data, c("treatment", "timepoint")))
})

test_that("validate_factors rejects non-dataframe input", {
  expect_error(
    validate_factors(list(a = 1), "a"),
    "must be a data frame"
  )
})

test_that("validate_factors rejects empty factors", {
  data <- data.frame(a = 1:3)
  expect_error(
    validate_factors(data, character(0)),
    "At least one factor"
  )
})


# Tests for validate_positive
test_that("validate_positive accepts positive values", {
  expect_true(validate_positive(c(1, 2, 3)))
  expect_true(validate_positive(c(0.001, 100, 1e6)))
})

test_that("validate_positive allows NA values", {
  expect_true(validate_positive(c(1, NA, 3)))
  expect_true(validate_positive(c(NA, NA, 1)))
})

test_that("validate_positive rejects zeros", {
  expect_error(
    validate_positive(c(1, 0, 3)),
    "Gamma GLM requires strictly positive values"
  )
  expect_error(
    validate_positive(c(0, 0, 0)),
    "3 zeros"
  )
})

test_that("validate_positive rejects negative values", {
  expect_error(
    validate_positive(c(1, -1, 3)),
    "1 negative"
  )
})

test_that("validate_positive includes variable name in error", {
  expect_error(
    validate_positive(c(0, 1, 2), "abundance"),
    "abundance"
  )
})


# Tests for check_zeros
test_that("check_zeros is an alias for validate_positive", {
  expect_true(check_zeros(c(1, 2, 3)))
  expect_error(check_zeros(c(0, 1, 2)), "strictly positive")
})


# Tests for compute_missingness
test_that("compute_missingness calculates NA rates correctly", {
  data <- tibble::tibble(
    peptide = rep(c("PEP_1", "PEP_2"), each = 4),
    value = c(1, 2, NA, 4, NA, NA, 3, 4)
  )

  result <- compute_missingness(data)

  expect_equal(nrow(result), 2)
  expect_equal(result$na_rate[result$peptide == "PEP_1"], 0.25)
  expect_equal(result$na_rate[result$peptide == "PEP_2"], 0.5)
})

test_that("compute_missingness calculates mean abundance", {
  data <- tibble::tibble(
    peptide = rep("PEP_1", 4),
    value = c(10, 20, 30, 40)
  )

  result <- compute_missingness(data)

  expect_equal(result$mean_abundance, 25)
})

test_that("compute_missingness includes MNAR score", {
  data <- tibble::tibble(
    peptide = rep(c("PEP_1", "PEP_2"), each = 4),
    value = c(1000, 1000, 1000, 1000, 10, NA, NA, NA)
  )

  result <- compute_missingness(data)

  expect_true("mnar_score" %in% names(result))
  # PEP_2 has more missing and lower abundance, should have higher MNAR score
  expect_true(result$mnar_score[result$peptide == "PEP_2"] >
                result$mnar_score[result$peptide == "PEP_1"])
})

test_that("compute_missingness requires peptide and value columns", {
  expect_error(
    compute_missingness(data.frame(a = 1)),
    "peptide.*value"
  )
})


# Tests for compute_design_summary
test_that("compute_design_summary counts replicates correctly", {
  data <- tibble::tibble(
    treatment = rep(c("ctrl", "trt"), each = 6),
    bio_rep = rep(1:3, 4),
    peptide = rep(c("PEP_1", "PEP_2"), 6),
    value = 1:12
  )

  result <- compute_design_summary(data, "treatment")

  expect_equal(nrow(result), 2)
  expect_true(all(result$n_reps == 3))
  expect_true(all(result$n_peptides == 2))
})

test_that("compute_design_summary handles multiple factors", {
  data <- make_factorial_test_data()

  result <- compute_design_summary(data, c("treatment", "timepoint"))

  expect_equal(nrow(result), 4)  # 2 treatments x 2 timepoints
  expect_true("treatment" %in% names(result))
  expect_true("timepoint" %in% names(result))
})


# Tests for mean_fold_change
test_that("mean_fold_change calculates correctly for vectors", {
  treatment <- c(20, 20, 20)
  control <- c(10, 10, 10)

  expect_equal(mean_fold_change(treatment, control), 2)
})

test_that("mean_fold_change handles NA values", {
  treatment <- c(20, NA, 20)
  control <- c(10, 10, NA)

  expect_equal(mean_fold_change(treatment, control), 2)
})

test_that("mean_fold_change works with matrices", {
  treatment <- matrix(c(20, 40, 20, 40), nrow = 2)
  control <- matrix(c(10, 20, 10, 20), nrow = 2)

  result <- mean_fold_change(treatment, control)

  expect_equal(result, c(2, 2))
})


# Tests for log2_fc
test_that("log2_fc calculates log2 fold change", {
  expect_equal(log2_fc(2), 1)
  expect_equal(log2_fc(4), 2)
  expect_equal(log2_fc(0.5), -1)
})


# Tests for build_formula
test_that("build_formula creates simple formula", {
  f <- build_formula("value", "treatment", interaction = FALSE)
  expect_equal(deparse(f), "value ~ treatment")
})

test_that("build_formula creates interaction formula", {
  f <- build_formula("value", c("treatment", "timepoint"), interaction = TRUE)
  expect_equal(deparse(f), "value ~ treatment * timepoint")
})

test_that("build_formula creates additive formula when no interaction", {
  f <- build_formula("value", c("treatment", "timepoint"), interaction = FALSE)
  expect_equal(deparse(f), "value ~ treatment + timepoint")
})


# Tests for is_useable
test_that("is_useable identifies NA and NaN", {
  x <- c(1, NA, NaN, 4)
  result <- is_useable(x)
  expect_equal(result, c(TRUE, FALSE, FALSE, TRUE))
})
