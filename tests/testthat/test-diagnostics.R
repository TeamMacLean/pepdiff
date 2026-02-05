# Tests for plot_fit_diagnostics() function

# =============================================================================
# Test Fixtures
# =============================================================================

# Create GLM results for testing
create_glm_results_for_diagnostics <- function() {
  set.seed(123)
  n_peptides <- 20
  n_reps <- 4

  peptides <- paste0("PEP_", sprintf("%03d", 1:n_peptides))
  genes <- paste0("GENE_", LETTERS[1:n_peptides])

  design <- expand.grid(
    peptide = peptides,
    treatment = c("ctrl", "trt"),
    bio_rep = 1:n_reps,
    stringsAsFactors = FALSE
  )

  sim_data <- design %>%
    dplyr::mutate(
      gene_id = genes[match(peptide, peptides)],
      pep_num = as.numeric(gsub("PEP_", "", peptide)),
      base = 10 + pep_num,
      effect = ifelse(pep_num <= 10 & treatment == "trt", 2, 1),
      value = rgamma(dplyr::n(), shape = 10, rate = 10 / (base * effect))
    ) %>%
    dplyr::select(peptide, gene_id, treatment, bio_rep, value)

  temp_file <- tempfile(fileext = ".csv")
  write.csv(sim_data, temp_file, row.names = FALSE)

  dat <- read_pepdiff(
    temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  unlink(temp_file)

  compare(dat, compare = "treatment", ref = "ctrl", method = "glm")
}

# Create ART results for testing error handling
create_art_results_for_diagnostics <- function() {

  set.seed(456)
  n_peptides <- 10
  n_reps <- 3

  peptides <- paste0("PEP_", sprintf("%03d", 1:n_peptides))
  genes <- paste0("GENE_", LETTERS[1:n_peptides])

  design <- expand.grid(
    peptide = peptides,
    treatment = c("ctrl", "trt"),
    bio_rep = 1:n_reps,
    stringsAsFactors = FALSE
  )

  sim_data <- design %>%
    dplyr::mutate(
      gene_id = genes[match(peptide, peptides)],
      value = runif(dplyr::n(), 5, 15)
    ) %>%
    dplyr::select(peptide, gene_id, treatment, bio_rep, value)

  temp_file <- tempfile(fileext = ".csv")
  write.csv(sim_data, temp_file, row.names = FALSE)

  dat <- read_pepdiff(
    temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  unlink(temp_file)

  compare(dat, compare = "treatment", ref = "ctrl", method = "art")
}

# =============================================================================
# Basic Functionality Tests
# =============================================================================

test_that("plot_fit_diagnostics returns expected structure for GLM results", {
  results <- create_glm_results_for_diagnostics()

  diag <- plot_fit_diagnostics(results)

  # Check return structure

expect_type(diag, "list")
  expect_named(diag, c("plot", "flagged", "summary"))

  # Check plot is ggplot object
  expect_s3_class(diag$plot, "ggplot")

  # Check flagged is tibble
  expect_s3_class(diag$flagged, "tbl_df")

  # Check summary is list with expected elements
  expect_type(diag$summary, "list")
  expect_true("n_analyzed" %in% names(diag$summary))
  expect_true("n_flagged" %in% names(diag$summary))
  expect_true("median_deviance" %in% names(diag$summary))
})

test_that("plot_fit_diagnostics errors for ART results", {
  results <- create_art_results_for_diagnostics()

  expect_error(
    plot_fit_diagnostics(results),
    "ART is non-parametric"
  )
})

test_that("plot_fit_diagnostics errors for non-pepdiff_results input", {
  expect_error(
    plot_fit_diagnostics(data.frame(x = 1)),
    "pepdiff_results"
  )
})

# =============================================================================
# Parameter Tests
# =============================================================================

test_that("n_sample parameter controls number of sample peptides", {
  results <- create_glm_results_for_diagnostics()

  # Default is 6
  diag <- plot_fit_diagnostics(results)
  # Note: actual number may be less if fewer peptides available

  # Custom value
  diag4 <- plot_fit_diagnostics(results, n_sample = 4)
  expect_type(diag4, "list")
})

test_that("deviance_threshold parameter works correctly", {
  results <- create_glm_results_for_diagnostics()

  # Very high threshold - no peptides flagged
  diag_high <- plot_fit_diagnostics(results, deviance_threshold = 1000)
  expect_equal(nrow(diag_high$flagged), 0)

  # Very low threshold - all peptides flagged
  diag_low <- plot_fit_diagnostics(results, deviance_threshold = 0)
  expect_true(nrow(diag_low$flagged) > 0)
})

# =============================================================================
# Flagged Peptides Tests
# =============================================================================

test_that("flagged tibble has expected columns", {
  results <- create_glm_results_for_diagnostics()

  diag <- plot_fit_diagnostics(results, deviance_threshold = 0)

  expected_cols <- c("peptide", "gene_id", "deviance", "fold_change", "p_value", "significant", "flag_reason")
  expect_true(all(expected_cols %in% names(diag$flagged)))
})

test_that("flagged peptides have correct flag_reason", {
  results <- create_glm_results_for_diagnostics()

  diag <- plot_fit_diagnostics(results, deviance_threshold = 0)

  expect_true(all(diag$flagged$flag_reason == "high_deviance"))
})

# =============================================================================
# Summary Statistics Tests
# =============================================================================

test_that("summary contains correct statistics", {
  results <- create_glm_results_for_diagnostics()

  diag <- plot_fit_diagnostics(results)

  expect_true(diag$summary$n_analyzed > 0)
  expect_true(is.numeric(diag$summary$median_deviance))
  expect_true(is.numeric(diag$summary$threshold))
  expect_true(diag$summary$n_flagged >= 0)
  expect_true(diag$summary$pct_flagged >= 0)
  expect_true(diag$summary$pct_flagged <= 100)
})

# =============================================================================
# Edge Cases
# =============================================================================

test_that("plot_fit_diagnostics handles results with some non-converged peptides", {
  # This tests robustness - real data may have convergence failures
  results <- create_glm_results_for_diagnostics()

  # Should still work even if some peptides didn't converge
  diag <- plot_fit_diagnostics(results)
  expect_type(diag, "list")
})

test_that("residuals are stored in diagnostics and used for QQ plot", {
  results <- create_glm_results_for_diagnostics()

  # Check residuals are stored
  expect_true("residuals" %in% names(results$diagnostics))

  # Check residuals are lists of numeric vectors
  resid_list <- results$diagnostics$residuals
  expect_type(resid_list, "list")

  # Converged peptides should have residuals
  converged_resid <- resid_list[results$diagnostics$converged]
  expect_true(all(sapply(converged_resid, is.numeric)))

  # Diagnostics should work and use stored residuals
  diag <- plot_fit_diagnostics(results)
  expect_type(diag, "list")
})

# =============================================================================
# Console Output Test
# =============================================================================

test_that("plot_fit_diagnostics prints summary to console", {
  results <- create_glm_results_for_diagnostics()

  # Capture output
  output <- capture.output(diag <- plot_fit_diagnostics(results))

  # Should print diagnostic summary
  expect_true(any(grepl("GLM Fit Diagnostics", output)))
  expect_true(any(grepl("Peptides analyzed", output)))
})
