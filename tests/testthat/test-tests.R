# Tests for tests.R - pairwise statistical tests

# =============================================================================
# Test Wilcoxon
# =============================================================================

test_that("test_wilcoxon returns p-value for different groups", {
  ctrl <- c(100, 110, 105, 108)
  trt <- c(200, 210, 205, 208)

  result <- test_wilcoxon(ctrl, trt)

  expect_type(result, "list")
  expect_true("p_value" %in% names(result))
  expect_true("method" %in% names(result))
  expect_equal(result$method, "wilcoxon")
  expect_true(result$p_value < 0.05)  # Groups are clearly different
})


test_that("test_wilcoxon returns high p-value for similar groups", {
  ctrl <- c(100, 110, 105, 108)
  trt <- c(102, 108, 107, 109)

  result <- test_wilcoxon(ctrl, trt)

  expect_true(result$p_value > 0.05)  # Groups are similar
})


test_that("test_wilcoxon handles NA values", {
  ctrl <- c(100, NA, 105, 108)
  trt <- c(200, 210, NA, 208)

  result <- test_wilcoxon(ctrl, trt)

  expect_false(is.na(result$p_value))
})


test_that("test_wilcoxon returns NA for insufficient data", {
  ctrl <- c(100)
  trt <- c(200, 210)

  result <- test_wilcoxon(ctrl, trt)

  expect_true(is.na(result$p_value))
})


# =============================================================================
# Test Bootstrap t-test
# =============================================================================

test_that("test_bootstrap_t returns p-value for different groups", {
  set.seed(123)
  ctrl <- c(100, 110, 105, 108, 103)
  trt <- c(200, 210, 205, 208, 203)

  result <- test_bootstrap_t(ctrl, trt, n_boot = 500, seed = 42)

  expect_type(result, "list")
  expect_true("p_value" %in% names(result))
  expect_true("t_obs" %in% names(result))
  expect_equal(result$method, "bootstrap_t")
  expect_true(result$p_value < 0.05)  # Groups are clearly different
})


test_that("test_bootstrap_t returns high p-value for similar groups", {
  ctrl <- c(100, 110, 105, 108, 103)
  trt <- c(102, 108, 107, 109, 104)

  result <- test_bootstrap_t(ctrl, trt, n_boot = 500, seed = 42)

  expect_true(result$p_value > 0.05)  # Groups are similar
})


test_that("test_bootstrap_t handles NA values", {
  ctrl <- c(100, NA, 105, 108, 103)
  trt <- c(200, 210, NA, 208, 203)

  result <- test_bootstrap_t(ctrl, trt, n_boot = 200, seed = 42)

  expect_false(is.na(result$p_value))
})


test_that("test_bootstrap_t is reproducible with seed", {
  ctrl <- c(100, 110, 105, 108, 103)
  trt <- c(150, 160, 155, 158, 153)

  result1 <- test_bootstrap_t(ctrl, trt, n_boot = 200, seed = 123)
  result2 <- test_bootstrap_t(ctrl, trt, n_boot = 200, seed = 123)

  expect_equal(result1$p_value, result2$p_value)
})


test_that("test_bootstrap_t returns NA for insufficient data", {
  ctrl <- c(100)
  trt <- c(200, 210)

  result <- test_bootstrap_t(ctrl, trt)

  expect_true(is.na(result$p_value))
})


# =============================================================================
# Test Bayes Factor t-test
# =============================================================================

test_that("test_bayes_t returns BF for different groups", {
  ctrl <- c(100, 110, 105, 108, 103)
  trt <- c(200, 210, 205, 208, 203)

  result <- test_bayes_t(ctrl, trt)

  expect_type(result, "list")
  expect_true("bf" %in% names(result))
  expect_true("effect_size" %in% names(result))
  expect_equal(result$method, "bayes_t")

  # BF should favor alternative for clearly different groups
  expect_true(result$bf > 1)
})


test_that("test_bayes_t returns small BF for similar groups", {
  ctrl <- c(100, 110, 105, 108, 103)
  trt <- c(102, 108, 107, 109, 104)

  result <- test_bayes_t(ctrl, trt)

  # BF should be close to 1 or favor null for similar groups
  expect_true(result$bf < 10)  # Not strong evidence for difference
})


test_that("test_bayes_t handles NA values", {
  ctrl <- c(100, NA, 105, 108, 103)
  trt <- c(200, 210, NA, 208, 203)

  result <- test_bayes_t(ctrl, trt)

  expect_false(is.na(result$bf))
  expect_false(is.na(result$effect_size))
})


test_that("test_bayes_t returns NA for insufficient data", {
  ctrl <- c(100)
  trt <- c(200, 210)

  result <- test_bayes_t(ctrl, trt)

  expect_true(is.na(result$bf))
})


test_that("test_bayes_t effect_size is positive when treatment > control", {
  ctrl <- c(100, 110, 105, 108, 103)
  trt <- c(200, 210, 205, 208, 203)

  result <- test_bayes_t(ctrl, trt)

  expect_true(result$effect_size > 0)
})


# =============================================================================
# Test Rank Products
# =============================================================================

test_that("test_rankprod returns p-values for different groups", {
  ctrl <- c(100, 110, 105, 108)
  trt <- c(200, 210, 205, 208)

  result <- test_rankprod(ctrl, trt, n_perm = 500, seed = 42)

  expect_type(result, "list")
  expect_true("p_value" %in% names(result))
  expect_true("p_value_up" %in% names(result))
  expect_true("p_value_down" %in% names(result))
  expect_equal(result$method, "rankprod")

  # p-values should be valid probabilities
  expect_true(result$p_value >= 0 && result$p_value <= 1)
  expect_true(result$p_value_up >= 0 && result$p_value_up <= 1)
  expect_true(result$p_value_down >= 0 && result$p_value_down <= 1)

  # Rank products should be positive
  expect_true(result$rp_up > 0)
  expect_true(result$rp_down > 0)
})


test_that("test_rankprod is reproducible with seed", {
  ctrl <- c(100, 110, 105, 108)
  trt <- c(150, 160, 155, 158)

  result1 <- test_rankprod(ctrl, trt, n_perm = 100, seed = 123)
  result2 <- test_rankprod(ctrl, trt, n_perm = 100, seed = 123)

  expect_equal(result1$p_value, result2$p_value)
})


test_that("test_rankprod handles NA values", {
  ctrl <- c(100, NA, 105, 108)
  trt <- c(200, 210, NA, 208)

  result <- test_rankprod(ctrl, trt, n_perm = 100, seed = 42)

  expect_false(is.na(result$p_value))
})


test_that("test_rankprod returns NA for insufficient data", {
  ctrl <- c(100)
  trt <- c(200, 210)

  result <- test_rankprod(ctrl, trt)

  expect_true(is.na(result$p_value))
})


# =============================================================================
# Test internal helper functions
# =============================================================================

test_that("calc_t_statistic computes correctly", {
  x <- c(10, 12, 11)
  y <- c(5, 6, 5.5)

  t_stat <- calc_t_statistic(x, y)

  # Should be positive since x > y
  expect_true(t_stat > 0)
})


test_that("calc_t_statistic returns NA for insufficient data", {
  x <- c(10)
  y <- c(5, 6)

  t_stat <- calc_t_statistic(x, y)

  expect_true(is.na(t_stat))
})


# =============================================================================
# Test legacy wrapper functions
# =============================================================================

test_that("get_bootstrap_percentile works with matrix input", {
  # Create matrix data (rows = peptides, cols = replicates)
  ctrl_matrix <- matrix(
    c(100, 110, 105, 200, 210, 205),
    nrow = 2, byrow = TRUE
  )
  trt_matrix <- matrix(
    c(150, 160, 155, 180, 190, 185),
    nrow = 2, byrow = TRUE
  )

  result <- get_bootstrap_percentile(trt_matrix, ctrl_matrix, iters = 100)

  expect_true("bootstrap_t_p_val" %in% names(result))
  expect_true("bootstrap_t_fdr" %in% names(result))
  expect_equal(nrow(result), 2)
})


test_that("get_wilcoxon_percentile works with matrix input", {
  ctrl_matrix <- matrix(
    c(100, 110, 105, 200, 210, 205),
    nrow = 2, byrow = TRUE
  )
  trt_matrix <- matrix(
    c(150, 160, 155, 180, 190, 185),
    nrow = 2, byrow = TRUE
  )

  result <- get_wilcoxon_percentile(trt_matrix, ctrl_matrix)

  expect_true("wilcoxon_p_val" %in% names(result))
  expect_true("wilcoxon_fdr" %in% names(result))
  expect_equal(nrow(result), 2)
})


test_that("get_kruskal_percentile works with matrix input", {
  ctrl_matrix <- matrix(
    c(100, 110, 105, 108, 200, 210, 205, 208),
    nrow = 2, byrow = TRUE
  )
  trt_matrix <- matrix(
    c(150, 160, 155, 158, 180, 190, 185, 188),
    nrow = 2, byrow = TRUE
  )

  result <- get_kruskal_percentile(trt_matrix, ctrl_matrix)

  expect_true("kruskal_p_val" %in% names(result))
  expect_true("kruskal_fdr" %in% names(result))
  expect_equal(nrow(result), 2)
})
