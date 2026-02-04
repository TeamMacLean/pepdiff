# Tests for Rank Products fix
# These tests verify that rankprod properly uses the RankProd package
# to rank across ALL peptides, not within single peptides.

# =============================================================================
# Test that rankprod detects true positives
# =============================================================================

test_that("rankprod detects truly differential peptides", {
  skip_if_not_installed("RankProd")

  # Create data with clear differential peptides
  # 30 peptides, 10 replicates per group, 10 truly differential (3-fold)
  set.seed(42)
  n_peptides <- 30
  n_reps <- 10

  peptides <- paste0("PEP_", sprintf("%03d", 1:n_peptides))
  genes <- paste0("GENE_", LETTERS[((1:n_peptides - 1) %% 26) + 1])

  # First 10 peptides are truly differential (5 up 3-fold, 5 down 0.33-fold)
  diff_peptides <- peptides[1:10]
  effect_directions <- c(rep(3, 5), rep(0.33, 5))

  sim_data <- expand.grid(
    peptide = peptides,
    treatment = c("ctrl", "trt"),
    bio_rep = 1:n_reps,
    stringsAsFactors = FALSE
  )
  sim_data$gene_id <- genes[match(sim_data$peptide, peptides)]
  sim_data$base <- rep(rgamma(n_peptides, shape = 5, rate = 0.5), each = 2 * n_reps)
  sim_data$effect <- ifelse(
    sim_data$peptide %in% diff_peptides & sim_data$treatment == "trt",
    effect_directions[match(sim_data$peptide, diff_peptides)],
    1
  )
  sim_data$value <- rgamma(nrow(sim_data), shape = 15, rate = 15 / (sim_data$base * sim_data$effect))
  sim_data <- sim_data[, c("peptide", "gene_id", "treatment", "bio_rep", "value")]

  temp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_file))
  readr::write_csv(sim_data, temp_file)

  dat <- read_pepdiff(
    temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  result <- compare(dat, compare = "treatment", ref = "ctrl",
                    method = "pairwise", test = "rankprod")

  # Should detect at least some of the 10 truly differential peptides
  sig_results <- result$results[result$results$significant == TRUE, ]
  sig_peptides <- sig_results$peptide

  # At least 5 of the 10 truly differential peptides should be detected
  n_true_positives <- sum(sig_peptides %in% diff_peptides)
  expect_gte(n_true_positives, 5,
             label = paste("true positives detected:", n_true_positives))
})


test_that("rankprod returns significant p-values for differential peptides",
{
  skip_if_not_installed("RankProd")

  # Simpler test: one clearly differential peptide
  set.seed(123)
  n_reps <- 5

  # Create data with one peptide that's 5-fold up
  sim_data <- data.frame(
    peptide = rep(c("PEP_DIFF", "PEP_NULL"), each = 2 * n_reps),
    gene_id = rep(c("GENE_A", "GENE_B"), each = 2 * n_reps),
    treatment = rep(c(rep("ctrl", n_reps), rep("trt", n_reps)), 2),
    bio_rep = rep(1:n_reps, 4),
    stringsAsFactors = FALSE
  )

  # PEP_DIFF: ctrl ~ 100, trt ~ 500 (5-fold up)
  # PEP_NULL: ctrl ~ 100, trt ~ 100 (no change)
  sim_data$value <- c(
    rnorm(n_reps, 100, 10),  # PEP_DIFF ctrl
    rnorm(n_reps, 500, 50),  # PEP_DIFF trt
    rnorm(n_reps, 100, 10),  # PEP_NULL ctrl
    rnorm(n_reps, 100, 10)   # PEP_NULL trt
  )

  temp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_file))
  readr::write_csv(sim_data, temp_file)

  dat <- read_pepdiff(
    temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  result <- compare(dat, compare = "treatment", ref = "ctrl",
                    method = "pairwise", test = "rankprod")

  # The differential peptide should have a low p-value
  diff_pval <- result$results$p_value[result$results$peptide == "PEP_DIFF"]
  null_pval <- result$results$p_value[result$results$peptide == "PEP_NULL"]

  expect_lt(diff_pval, 0.1,
            label = paste("differential peptide p-value:", diff_pval))

  # The null peptide should have a higher p-value
  expect_gt(null_pval, diff_pval,
            label = "null peptide p-value vs differential")
})


# =============================================================================
# Test missing RankProd package handling
# =============================================================================

test_that("rankprod gives informative warning when RankProd not installed", {
  # We can't actually test this without uninstalling RankProd

  # But we can test the message content by mocking
  skip("Cannot test missing package without uninstalling - manual verification needed")
})


test_that("rankprod returns NA p-values when RankProd not installed", {
  # Similar skip - need mock or manual test

  skip("Cannot test missing package without uninstalling - manual verification needed")
})


# =============================================================================
# Test results format
# =============================================================================

test_that("rankprod results have same structure as other pairwise tests", {
  skip_if_not_installed("RankProd")

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

  # Standard columns should be present
  expect_true("peptide" %in% names(result$results))
  expect_true("gene_id" %in% names(result$results))
  expect_true("comparison" %in% names(result$results))
  expect_true("fold_change" %in% names(result$results))
  expect_true("log2_fc" %in% names(result$results))
  expect_true("test" %in% names(result$results))
  expect_true("p_value" %in% names(result$results))
  expect_true("fdr" %in% names(result$results))
  expect_true("significant" %in% names(result$results))

  # test column should say "rankprod"
  expect_equal(unique(result$results$test), "rankprod")
})


test_that("rankprod FDR correction is applied", {
  skip_if_not_installed("RankProd")

  test_data <- make_test_data(n_peptides = 20, n_reps = 4, seed = 42)
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

  # FDR values should be >= p-values (where both are not NA)
  valid_idx <- !is.na(result$results$p_value) & !is.na(result$results$fdr)
  if (any(valid_idx)) {
    expect_true(all(result$results$fdr[valid_idx] >= result$results$p_value[valid_idx]))
  }
})


# =============================================================================
# Test integration with results methods
# =============================================================================

test_that("significant() works with rankprod results", {
  skip_if_not_installed("RankProd")

  # Use data with known differential peptides
  set.seed(42)
  n_peptides <- 20
  n_reps <- 6

  peptides <- paste0("PEP_", sprintf("%03d", 1:n_peptides))
  genes <- paste0("GENE_", LETTERS[((1:n_peptides - 1) %% 26) + 1])

  # First 5 peptides are 3-fold up
  diff_peptides <- peptides[1:5]

  sim_data <- expand.grid(
    peptide = peptides,
    treatment = c("ctrl", "trt"),
    bio_rep = 1:n_reps,
    stringsAsFactors = FALSE
  )
  sim_data$gene_id <- genes[match(sim_data$peptide, peptides)]
  sim_data$base <- rep(rgamma(n_peptides, shape = 5, rate = 0.5), each = 2 * n_reps)
  sim_data$effect <- ifelse(
    sim_data$peptide %in% diff_peptides & sim_data$treatment == "trt",
    3,  # 3-fold up
    1
  )
  sim_data$value <- rgamma(nrow(sim_data), shape = 15, rate = 15 / (sim_data$base * sim_data$effect))
  sim_data <- sim_data[, c("peptide", "gene_id", "treatment", "bio_rep", "value")]

  temp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_file))
  readr::write_csv(sim_data, temp_file)

  dat <- read_pepdiff(
    temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  result <- compare(dat, compare = "treatment", ref = "ctrl",
                    method = "pairwise", test = "rankprod")

  sig <- significant(result)

  # Should return a tibble

  expect_s3_class(sig, "tbl_df")

  # All returned results should have fdr < 0.05
  if (nrow(sig) > 0) {
    expect_true(all(sig$fdr < 0.05, na.rm = TRUE))
  }
})


test_that("print and summary work with rankprod results", {
  skip_if_not_installed("RankProd")

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

  # print should work without error
  output <- capture.output(print(result))
  expect_true(length(output) > 0)

  # summary should work without error
  summ_output <- capture.output(summ <- summary(result))
  expect_true(length(summ_output) > 0)
  expect_type(summ, "list")
})


test_that("plot works with rankprod results", {
  skip_if_not_installed("RankProd")

  test_data <- make_test_data(n_peptides = 10, n_reps = 4, seed = 42)
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

  # plot should return a ggplot/grob without error
  p <- plot(result)
  expect_true(!is.null(p))
})


# =============================================================================
# Test deprecation of test_rankprod()
# =============================================================================

test_that("test_rankprod() shows deprecation warning", {
  # This test expects the deprecated function to warn users
  # It should pass after we add .Deprecated() call
  expect_warning(
    test_rankprod(c(1, 2, 3, 4), c(5, 6, 7, 8), n_perm = 10),
    "deprecated|Deprecated"
  )
})


# =============================================================================
# Test that rankprod is handled differently from other tests
# =============================================================================

test_that("rankprod is dispatched separately from per-peptide tests", {
  skip_if_not_installed("RankProd")

  # This test verifies that rankprod uses the full matrix approach
  # by checking that it produces different (better) results than
  # the broken per-peptide approach would

  set.seed(42)
  n_peptides <- 10
  n_reps <- 5

  peptides <- paste0("PEP_", sprintf("%03d", 1:n_peptides))
  genes <- paste0("GENE_", LETTERS[1:n_peptides])

  # Create data where peptide 1 is clearly up-regulated
  sim_data <- expand.grid(
    peptide = peptides,
    treatment = c("ctrl", "trt"),
    bio_rep = 1:n_reps,
    stringsAsFactors = FALSE
  )
  sim_data$gene_id <- genes[match(sim_data$peptide, peptides)]

  # Peptide 1: 3-fold up; others: no change
  sim_data$value <- mapply(function(pep, trt) {
    if (pep == "PEP_001" && trt == "trt") {
      rnorm(1, 300, 30)
    } else {
      rnorm(1, 100, 15)
    }
  }, sim_data$peptide, sim_data$treatment)

  temp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_file))
  readr::write_csv(sim_data, temp_file)

  dat <- read_pepdiff(
    temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  result <- compare(dat, compare = "treatment", ref = "ctrl",
                    method = "pairwise", test = "rankprod")

  # PEP_001 should have the smallest p-value
  pep_001_pval <- result$results$p_value[result$results$peptide == "PEP_001"]
  other_pvals <- result$results$p_value[result$results$peptide != "PEP_001"]

  expect_lt(pep_001_pval, min(other_pvals, na.rm = TRUE),
            label = paste("differential peptide p-value:", pep_001_pval, "vs min other:", min(other_pvals, na.rm = TRUE)))

  # P-value should be reasonably small (not ~1 as broken implementation gives)
  expect_lt(pep_001_pval, 0.5,
            label = paste("p-value:", pep_001_pval))
})
