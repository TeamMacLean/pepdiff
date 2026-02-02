# pepdiff_results S3 class for analysis results
# Contains differential abundance results in long format

# =============================================================================
# pepdiff_results Class Constructor and Validator
# =============================================================================

#' Create a new pepdiff_results object
#'
#' Low-level constructor for pepdiff_results objects. Typically created by
#' [compare()] rather than directly.
#'
#' @param results Tibble with results in long format
#' @param comparisons Tibble defining the comparisons made
#' @param method Character, the statistical method used
#' @param diagnostics Tibble with model diagnostics (convergence, etc.)
#' @param params List of parameters used in the analysis
#' @param data The original pepdiff_data object
#' @param call The original function call
#'
#' @return A pepdiff_results object
#' @keywords internal
new_pepdiff_results <- function(results, comparisons, method, diagnostics, params, data, call) {
  structure(
    list(
      results = results,
      comparisons = comparisons,
      method = method,
      diagnostics = diagnostics,
      params = params,
      data = data,
      call = call
    ),
    class = "pepdiff_results"
  )
}


#' Validate a pepdiff_results object
#'
#' @param x A pepdiff_results object to validate
#' @return The validated object (invisibly), or throws an error
#' @keywords internal
validate_pepdiff_results <- function(x) {
  if (!inherits(x, "pepdiff_results")) {
    stop("Object must be of class 'pepdiff_results'", call. = FALSE)
  }

  # Check required components
  required <- c("results", "comparisons", "method", "diagnostics", "params", "data", "call")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(
      "pepdiff_results missing required components: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Check results is a tibble/data.frame
  if (!is.data.frame(x$results)) {
    stop("'results' component must be a data frame", call. = FALSE)
  }

  # Check required columns in results
  required_cols <- c("peptide", "gene_id", "comparison", "p_value")
  missing_cols <- setdiff(required_cols, names(x$results))
  if (length(missing_cols) > 0) {
    stop(
      "Results missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Check method is valid
  valid_methods <- c("glm", "art", "pairwise")
  if (!x$method %in% valid_methods) {
    stop(
      "Method must be one of: ",
      paste(valid_methods, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(x)
}


# =============================================================================
# Print Method
# =============================================================================

#' Print method for pepdiff_results
#'
#' @param x A pepdiff_results object
#' @param ... Additional arguments (ignored)
#' @return The object invisibly
#' @export
print.pepdiff_results <- function(x, ...) {
  cat("pepdiff_results object\n")
  cat("----------------------\n")

  # Method
  cat(sprintf("Method: %s\n", x$method))

  # Summary counts
  n_peptides <- length(unique(x$results$peptide))
  n_comparisons <- length(unique(x$results$comparison))
  n_tests <- nrow(x$results)

  cat(sprintf("Peptides: %d\n", n_peptides))
  cat(sprintf("Comparisons: %d\n", n_comparisons))
  cat(sprintf("Total tests: %d\n", n_tests))

  # Significance summary
  alpha <- x$params$alpha %||% 0.05
  if ("fdr" %in% names(x$results)) {
    n_sig <- sum(x$results$fdr < alpha, na.rm = TRUE)
    cat(sprintf("\nSignificant (FDR < %.2f): %d (%.1f%%)\n",
                alpha, n_sig, 100 * n_sig / n_tests))
  }

  if ("significant" %in% names(x$results)) {
    n_sig <- sum(x$results$significant, na.rm = TRUE)
    cat(sprintf("Marked significant: %d\n", n_sig))
  }

  # Convergence warnings
  if (!is.null(x$diagnostics) && "converged" %in% names(x$diagnostics)) {
    n_failed <- sum(!x$diagnostics$converged, na.rm = TRUE)
    if (n_failed > 0) {
      cat(sprintf("\nWarning: %d peptides excluded (model did not converge)\n", n_failed))
    }
  }

  invisible(x)
}


# =============================================================================
# Summary Method
# =============================================================================

#' Summary method for pepdiff_results
#'
#' @param object A pepdiff_results object
#' @param ... Additional arguments (ignored)
#' @return A summary list invisibly
#' @export
summary.pepdiff_results <- function(object, ...) {
  cat("pepdiff_results Summary\n")
  cat("=======================\n\n")

  # Method and parameters
  cat("Analysis method:", object$method, "\n")
  if (!is.null(object$params$test)) {
    cat("Test:", object$params$test, "\n")
  }
  cat("Alpha:", object$params$alpha %||% 0.05, "\n")
  cat("FDR method:", object$params$fdr_method %||% "BH", "\n\n")

  # Per-comparison breakdown
  cat("Results by comparison:\n")
  alpha <- object$params$alpha %||% 0.05

  comparison_summary <- object$results %>%
    dplyr::group_by(.data$comparison) %>%
    dplyr::summarize(
      n_peptides = dplyr::n(),
      n_tested = sum(!is.na(.data$p_value)),
      n_significant = if ("fdr" %in% names(object$results)) sum(.data$fdr < alpha, na.rm = TRUE) else NA_integer_,
      pct_significant = if ("fdr" %in% names(object$results)) 100 * mean(.data$fdr < alpha, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )

  print(comparison_summary, n = Inf)

  # Convergence summary
  if (!is.null(object$diagnostics) && "converged" %in% names(object$diagnostics)) {
    n_total <- nrow(object$diagnostics)
    n_converged <- sum(object$diagnostics$converged, na.rm = TRUE)
    n_failed <- n_total - n_converged

    cat("\nModel convergence:\n")
    cat(sprintf("  Converged: %d (%.1f%%)\n", n_converged, 100 * n_converged / n_total))
    if (n_failed > 0) {
      cat(sprintf("  Failed: %d (%.1f%%)\n", n_failed, 100 * n_failed / n_total))
    }
  }

  invisible(list(
    method = object$method,
    params = object$params,
    comparison_summary = comparison_summary,
    n_peptides = length(unique(object$results$peptide)),
    n_comparisons = length(unique(object$results$comparison))
  ))
}


# =============================================================================
# Accessor Functions
# =============================================================================

#' Extract significant results
#'
#' @param x A pepdiff_results object
#' @param alpha Significance threshold (default uses analysis alpha)
#' @param by_fdr Logical, use FDR-adjusted p-values (default TRUE)
#'
#' @return A tibble of significant results
#' @export
significant <- function(x, alpha = NULL, by_fdr = TRUE) {
  UseMethod("significant")
}


#' @export
significant.pepdiff_results <- function(x, alpha = NULL, by_fdr = TRUE) {
  if (is.null(alpha)) {
    alpha <- x$params$alpha %||% 0.05
  }

  if (by_fdr && "fdr" %in% names(x$results)) {
    x$results[x$results$fdr < alpha & !is.na(x$results$fdr), ]
  } else if ("p_value" %in% names(x$results)) {
    x$results[x$results$p_value < alpha & !is.na(x$results$p_value), ]
  } else {
    x$results[0, ]  # Empty tibble with same structure
  }
}


#' Get results for a specific peptide
#'
#' @param x A pepdiff_results object
#' @param peptide Peptide ID to retrieve
#'
#' @return A tibble with results for the specified peptide
#' @export
get_peptide <- function(x, peptide) {
  UseMethod("get_peptide")
}


#' @export
get_peptide.pepdiff_results <- function(x, peptide) {
  x$results[x$results$peptide == peptide, ]
}


#' Get results for a specific comparison
#'
#' @param x A pepdiff_results object
#' @param comparison Comparison name to retrieve
#'
#' @return A tibble with results for the specified comparison
#' @export
get_comparison <- function(x, comparison) {
  UseMethod("get_comparison")
}


#' @export
get_comparison.pepdiff_results <- function(x, comparison) {
  x$results[x$results$comparison == comparison, ]
}


# =============================================================================
# Results Building Helper
# =============================================================================

#' Build results tibble from model output
#'
#' Combines model results into the standard long-format results tibble.
#'
#' @param model_results List of results from run_models()
#' @param data Original pepdiff_data object
#' @param compare Factor being compared
#' @param method Method used
#' @param alpha Significance threshold
#' @param fdr_method FDR correction method
#'
#' @return A tibble in long format
#' @keywords internal
#' @importFrom rlang .data
build_results_tibble <- function(model_results, data, compare, method, alpha, fdr_method) {
  # Extract contrasts from each peptide
  results_list <- lapply(names(model_results), function(pep) {
    res <- model_results[[pep]]

    if (!res$converged || is.null(res$contrasts)) {
      return(NULL)
    }

    contrasts <- res$contrasts

    # Get gene_id for this peptide
    gene_id <- unique(data$data$gene_id[data$data$peptide == pep])[1]

    tibble::tibble(
      peptide = pep,
      gene_id = gene_id,
      comparison = contrasts$contrast,
      fold_change = if ("fold_change" %in% names(contrasts)) contrasts$fold_change else NA_real_,
      log2_fc = if ("fold_change" %in% names(contrasts)) log2(contrasts$fold_change) else NA_real_,
      estimate = if ("estimate" %in% names(contrasts)) contrasts$estimate else NA_real_,
      se = if ("se" %in% names(contrasts)) contrasts$se else NA_real_,
      test = method,
      p_value = contrasts$p_value
    )
  })

  # Combine all results
  results <- dplyr::bind_rows(results_list)

  if (nrow(results) == 0) {
    # Return empty tibble with correct structure
    return(tibble::tibble(
      peptide = character(),
      gene_id = character(),
      comparison = character(),
      fold_change = numeric(),
      log2_fc = numeric(),
      estimate = numeric(),
      se = numeric(),
      test = character(),
      p_value = numeric(),
      fdr = numeric(),
      significant = logical()
    ))
  }

  # Apply FDR correction within each comparison
  results <- results %>%
    dplyr::group_by(.data$comparison) %>%
    dplyr::mutate(
      fdr = stats::p.adjust(.data$p_value, method = fdr_method)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      significant = .data$fdr < alpha
    )

  results
}


#' Build results tibble from pairwise tests
#'
#' @param test_results List of test results for each peptide
#' @param data Original pepdiff_data object
#' @param compare Factor being compared
#' @param ref Reference level
#' @param test Test method name
#' @param alpha Significance threshold
#' @param fdr_method FDR correction method
#'
#' @return A tibble in long format
#' @keywords internal
#' @importFrom rlang .data
build_pairwise_results <- function(test_results, data, compare, ref, test, alpha, fdr_method) {
  # Get levels of the comparison factor
  levels <- unique(data$data[[compare]])
  treatment_levels <- setdiff(levels, ref)

  # Build comparison name
  comparison_name <- paste(treatment_levels[1], "vs", ref)

  results <- tibble::tibble(
    peptide = names(test_results),
    gene_id = vapply(names(test_results), function(pep) {
      unique(data$data$gene_id[data$data$peptide == pep])[1]
    }, character(1)),
    comparison = comparison_name,
    fold_change = vapply(test_results, function(x) x$fold_change %||% NA_real_, numeric(1)),
    log2_fc = vapply(test_results, function(x) {
      fc <- x$fold_change %||% NA_real_
      if (is.na(fc)) NA_real_ else log2(fc)
    }, numeric(1)),
    test = test,
    p_value = vapply(test_results, function(x) x$p_value %||% NA_real_, numeric(1))
  )

  # Apply FDR correction
  results <- results %>%
    dplyr::mutate(
      fdr = stats::p.adjust(.data$p_value, method = fdr_method),
      significant = .data$fdr < alpha
    )

  results
}


# =============================================================================
# Null coalescing operator
# =============================================================================

# Use rlang's null coalescing operator internally
# Note: %||% is imported from rlang
