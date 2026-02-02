# Tests for results.R - pepdiff_results class

# =============================================================================
# Test Constructor and Validator
# =============================================================================

test_that("new_pepdiff_results creates object with correct structure", {
  results <- tibble::tibble(
    peptide = c("PEP_001", "PEP_002"),
    gene_id = c("GENE_001", "GENE_001"),
    comparison = c("trt vs ctrl", "trt vs ctrl"),
    fold_change = c(2.0, 1.5),
    log2_fc = c(1.0, 0.585),
    test = c("glm", "glm"),
    p_value = c(0.01, 0.1),
    fdr = c(0.02, 0.1),
    significant = c(TRUE, FALSE)
  )

  comparisons <- tibble::tibble(
    comparison = "trt vs ctrl",
    contrast = "treatment",
    ref = "ctrl"
  )

  diagnostics <- tibble::tibble(
    peptide = c("PEP_001", "PEP_002"),
    converged = c(TRUE, TRUE)
  )

  params <- list(alpha = 0.05, fdr_method = "BH")

  # Create minimal mock pepdiff_data
  mock_data <- structure(
    list(
      data = tibble::tibble(
        peptide = c("PEP_001", "PEP_002"),
        gene_id = c("GENE_001", "GENE_001"),
        treatment = c("ctrl", "trt"),
        bio_rep = c("1", "1"),
        value = c(100, 200)
      ),
      factors = "treatment",
      peptides = c("PEP_001", "PEP_002")
    ),
    class = "pepdiff_data"
  )

  obj <- new_pepdiff_results(
    results = results,
    comparisons = comparisons,
    method = "glm",
    diagnostics = diagnostics,
    params = params,
    data = mock_data,
    call = quote(compare(data))
  )

  expect_s3_class(obj, "pepdiff_results")
  expect_equal(obj$method, "glm")
})


test_that("validate_pepdiff_results catches missing components", {
  bad_obj <- structure(
    list(results = tibble::tibble(), method = "glm"),
    class = "pepdiff_results"
  )

  expect_error(validate_pepdiff_results(bad_obj), "missing required components")
})


test_that("validate_pepdiff_results catches invalid method", {
  results <- tibble::tibble(
    peptide = "PEP_001",
    gene_id = "GENE_001",
    comparison = "trt vs ctrl",
    p_value = 0.05
  )

  bad_obj <- structure(
    list(
      results = results,
      comparisons = tibble::tibble(),
      method = "invalid_method",  # Invalid!
      diagnostics = tibble::tibble(),
      params = list(),
      data = NULL,
      call = quote(test())
    ),
    class = "pepdiff_results"
  )

  expect_error(validate_pepdiff_results(bad_obj), "Method must be one of")
})


# =============================================================================
# Test Print Method
# =============================================================================

test_that("print.pepdiff_results produces output", {
  results <- tibble::tibble(
    peptide = c("PEP_001", "PEP_002", "PEP_003"),
    gene_id = c("GENE_001", "GENE_001", "GENE_002"),
    comparison = rep("trt vs ctrl", 3),
    fold_change = c(2.0, 1.5, 0.5),
    log2_fc = c(1.0, 0.585, -1.0),
    test = rep("glm", 3),
    p_value = c(0.01, 0.1, 0.001),
    fdr = c(0.02, 0.15, 0.003),
    significant = c(TRUE, FALSE, TRUE)
  )

  obj <- new_pepdiff_results(
    results = results,
    comparisons = tibble::tibble(comparison = "trt vs ctrl"),
    method = "glm",
    diagnostics = tibble::tibble(peptide = c("PEP_001", "PEP_002", "PEP_003"), converged = c(TRUE, TRUE, TRUE)),
    params = list(alpha = 0.05, fdr_method = "BH"),
    data = NULL,
    call = quote(compare(data))
  )

  output <- capture.output(print(obj))

  expect_true(any(grepl("pepdiff_results", output)))
  expect_true(any(grepl("Method: glm", output)))
  expect_true(any(grepl("Significant", output)))
})


test_that("print shows convergence warnings", {
  results <- tibble::tibble(
    peptide = c("PEP_001", "PEP_002"),
    gene_id = c("GENE_001", "GENE_001"),
    comparison = rep("trt vs ctrl", 2),
    p_value = c(0.01, 0.1)
  )

  obj <- new_pepdiff_results(
    results = results,
    comparisons = tibble::tibble(comparison = "trt vs ctrl"),
    method = "glm",
    diagnostics = tibble::tibble(peptide = c("PEP_001", "PEP_002", "PEP_003"), converged = c(TRUE, TRUE, FALSE)),
    params = list(alpha = 0.05),
    data = NULL,
    call = quote(test())
  )

  output <- capture.output(print(obj))

  expect_true(any(grepl("Warning", output)))
  expect_true(any(grepl("did not converge", output)))
})


# =============================================================================
# Test Summary Method
# =============================================================================

test_that("summary.pepdiff_results produces detailed output", {
  results <- tibble::tibble(
    peptide = c("PEP_001", "PEP_002", "PEP_003", "PEP_004"),
    gene_id = rep("GENE_001", 4),
    comparison = c(rep("comp1", 2), rep("comp2", 2)),
    p_value = c(0.01, 0.1, 0.001, 0.05),
    fdr = c(0.02, 0.15, 0.002, 0.08)
  )

  obj <- new_pepdiff_results(
    results = results,
    comparisons = tibble::tibble(comparison = c("comp1", "comp2")),
    method = "glm",
    diagnostics = tibble::tibble(peptide = paste0("PEP_00", 1:4), converged = rep(TRUE, 4)),
    params = list(alpha = 0.05, fdr_method = "BH"),
    data = NULL,
    call = quote(test())
  )

  output <- capture.output(result <- summary(obj))

  expect_true(any(grepl("Summary", output)))
  expect_true(any(grepl("comparison", output)))
  expect_type(result, "list")
})


# =============================================================================
# Test Accessor Functions
# =============================================================================

test_that("significant extracts significant results", {
  results <- tibble::tibble(
    peptide = c("PEP_001", "PEP_002", "PEP_003"),
    gene_id = rep("GENE_001", 3),
    comparison = rep("trt vs ctrl", 3),
    p_value = c(0.01, 0.1, 0.001),
    fdr = c(0.02, 0.15, 0.003)
  )

  obj <- new_pepdiff_results(
    results = results,
    comparisons = tibble::tibble(comparison = "trt vs ctrl"),
    method = "glm",
    diagnostics = tibble::tibble(),
    params = list(alpha = 0.05),
    data = NULL,
    call = quote(test())
  )

  sig <- significant(obj)

  expect_equal(nrow(sig), 2)  # PEP_001 and PEP_003
  expect_true(all(sig$fdr < 0.05))
})


test_that("significant respects custom alpha", {
  results <- tibble::tibble(
    peptide = c("PEP_001", "PEP_002"),
    gene_id = rep("GENE_001", 2),
    comparison = rep("trt vs ctrl", 2),
    p_value = c(0.01, 0.1),
    fdr = c(0.02, 0.15)
  )

  obj <- new_pepdiff_results(
    results = results,
    comparisons = tibble::tibble(comparison = "trt vs ctrl"),
    method = "glm",
    diagnostics = tibble::tibble(),
    params = list(alpha = 0.05),
    data = NULL,
    call = quote(test())
  )

  # With strict alpha
  sig_strict <- significant(obj, alpha = 0.01)
  expect_equal(nrow(sig_strict), 0)

  # With lenient alpha
  sig_lenient <- significant(obj, alpha = 0.20)
  expect_equal(nrow(sig_lenient), 2)
})


test_that("get_peptide extracts single peptide results", {
  results <- tibble::tibble(
    peptide = c("PEP_001", "PEP_002", "PEP_001"),
    gene_id = rep("GENE_001", 3),
    comparison = c("comp1", "comp1", "comp2"),
    p_value = c(0.01, 0.1, 0.02)
  )

  obj <- new_pepdiff_results(
    results = results,
    comparisons = tibble::tibble(comparison = c("comp1", "comp2")),
    method = "glm",
    diagnostics = tibble::tibble(),
    params = list(),
    data = NULL,
    call = quote(test())
  )

  pep_results <- get_peptide(obj, "PEP_001")

  expect_equal(nrow(pep_results), 2)  # Two comparisons for PEP_001
  expect_true(all(pep_results$peptide == "PEP_001"))
})


test_that("get_comparison extracts single comparison results", {
  results <- tibble::tibble(
    peptide = c("PEP_001", "PEP_002", "PEP_001"),
    gene_id = rep("GENE_001", 3),
    comparison = c("comp1", "comp1", "comp2"),
    p_value = c(0.01, 0.1, 0.02)
  )

  obj <- new_pepdiff_results(
    results = results,
    comparisons = tibble::tibble(comparison = c("comp1", "comp2")),
    method = "glm",
    diagnostics = tibble::tibble(),
    params = list(),
    data = NULL,
    call = quote(test())
  )

  comp_results <- get_comparison(obj, "comp1")

  expect_equal(nrow(comp_results), 2)  # Two peptides for comp1
  expect_true(all(comp_results$comparison == "comp1"))
})


# =============================================================================
# Test Results Building Helpers
# =============================================================================

test_that("build_pairwise_results creates correct structure", {
  # Create mock test results
  test_results <- list(
    "PEP_001" = list(p_value = 0.01, fold_change = 2.0),
    "PEP_002" = list(p_value = 0.1, fold_change = 1.2)
  )

  # Create mock pepdiff_data
  mock_data <- structure(
    list(
      data = tibble::tibble(
        peptide = c("PEP_001", "PEP_001", "PEP_002", "PEP_002"),
        gene_id = c("GENE_001", "GENE_001", "GENE_001", "GENE_001"),
        treatment = c("ctrl", "trt", "ctrl", "trt"),
        bio_rep = c("1", "1", "1", "1"),
        value = c(100, 200, 100, 120)
      ),
      factors = "treatment",
      peptides = c("PEP_001", "PEP_002")
    ),
    class = "pepdiff_data"
  )

  results <- build_pairwise_results(
    test_results = test_results,
    data = mock_data,
    compare = "treatment",
    ref = "ctrl",
    test = "wilcoxon",
    alpha = 0.05,
    fdr_method = "BH"
  )

  expect_s3_class(results, "tbl_df")
  expect_true("peptide" %in% names(results))
  expect_true("gene_id" %in% names(results))
  expect_true("comparison" %in% names(results))
  expect_true("fold_change" %in% names(results))
  expect_true("p_value" %in% names(results))
  expect_true("fdr" %in% names(results))
  expect_true("significant" %in% names(results))
})


test_that("build_results_tibble handles empty input", {
  mock_data <- structure(
    list(
      data = tibble::tibble(
        peptide = character(),
        gene_id = character(),
        treatment = character(),
        bio_rep = character(),
        value = numeric()
      ),
      factors = "treatment",
      peptides = character()
    ),
    class = "pepdiff_data"
  )

  results <- build_results_tibble(
    model_results = list(),
    data = mock_data,
    compare = "treatment",
    method = "glm",
    alpha = 0.05,
    fdr_method = "BH"
  )

  expect_s3_class(results, "tbl_df")
  expect_equal(nrow(results), 0)
})


# =============================================================================
# Test %||% operator
# =============================================================================

test_that("null coalescing operator works correctly", {
  expect_equal(NULL %||% "default", "default")
  expect_equal("value" %||% "default", "value")
  expect_equal(5 %||% 10, 5)
})
