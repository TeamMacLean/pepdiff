# Main compare() function for differential abundance analysis
# This is the primary user-facing analysis function

# =============================================================================
# Generic and Methods
# =============================================================================

#' Compare peptide abundances between conditions
#'
#' Performs differential abundance analysis on proteomics data. Supports three
#' methods: GLM (default), ART (non-parametric), and pairwise tests.
#'
#' @param data A pepdiff_data object from [read_pepdiff()]
#' @param compare Factor to compare (character string)
#' @param ref Reference level for comparisons
#' @param within Optional factor(s) to stratify by
#' @param method Analysis method: "glm" (default), "art", or "pairwise"
#' @param test For pairwise method: "wilcoxon", "bootstrap_t", "bayes_t", or "rankprod"
#' @param alpha Significance threshold (default 0.05)
#' @param fdr_method FDR correction method (default "BH")
#' @param ... Additional arguments passed to methods
#'
#' @return A pepdiff_results object containing:
#'   \item{results}{Tibble with peptide, gene_id, comparison, fold_change, log2_fc, p_value, fdr, significant}
#'   \item{comparisons}{Tibble defining the comparisons made}
#'   \item{method}{Statistical method used}
#'   \item{diagnostics}{Model convergence information (for GLM/ART)}
#'   \item{params}{Analysis parameters}
#'   \item{data}{The original pepdiff_data object}
#'   \item{call}{The function call}
#'
#' @export
#' @examples
#' \dontrun{
#' # Simple comparison
#' results <- compare(data, compare = "treatment", ref = "ctrl")
#'
#' # Stratified comparison
#' results <- compare(data, compare = "treatment", ref = "ctrl", within = "timepoint")
#'
#' # Pairwise test
#' results <- compare(data, compare = "treatment", ref = "ctrl",
#'                    method = "pairwise", test = "wilcoxon")
#' }
compare <- function(data, ...) {
  UseMethod("compare")
}


#' @rdname compare
#' @export
#' @method compare pepdiff_data
compare.pepdiff_data <- function(data, compare, ref,
                                  within = NULL,
                                  method = c("glm", "art", "pairwise"),
                                  test = c("wilcoxon", "bootstrap_t", "bayes_t", "rankprod"),
                                  alpha = 0.05,
                                  fdr_method = "BH",
                                  ...) {
  # Capture call
  call <- match.call()

  # Validate inputs
  method <- match.arg(method)
  test <- match.arg(test)

  if (missing(compare) || !is.character(compare)) {
    stop("'compare' must be a character string naming the factor to compare", call. = FALSE)
  }

  if (!compare %in% data$factors) {
    stop(
      "Factor '", compare, "' not found. Available factors: ",
      paste(data$factors, collapse = ", "),
      call. = FALSE
    )
  }

  if (missing(ref)) {
    # Use first level as reference
    ref <- sort(unique(data$data[[compare]]))[1]
    message("Using '", ref, "' as reference level")
  }

  if (!ref %in% data$data[[compare]]) {
    stop("Reference level '", ref, "' not found in factor '", compare, "'", call. = FALSE)
  }

  if (!is.null(within) && !all(within %in% data$factors)) {
    missing_within <- setdiff(within, data$factors)
    stop(
      "'within' factor(s) not found: ",
      paste(missing_within, collapse = ", "),
      call. = FALSE
    )
  }

  # Dispatch to appropriate method
  if (method == "pairwise") {
    results <- compare_pairwise(data, compare, ref, within, test, alpha, fdr_method)
  } else if (method == "glm") {
    results <- compare_glm(data, compare, ref, within, alpha, fdr_method)
  } else if (method == "art") {
    results <- compare_art(data, compare, ref, within, alpha, fdr_method)
  }

  # Build comparisons tibble
  levels <- setdiff(unique(data$data[[compare]]), ref)
  comparisons <- tibble::tibble(
    comparison = paste(levels, "vs", ref),
    factor = compare,
    treatment = levels,
    reference = ref
  )

  # Construct pepdiff_results object
  result_obj <- new_pepdiff_results(
    results = results$results,
    comparisons = comparisons,
    method = method,
    diagnostics = results$diagnostics,
    params = list(
      compare = compare,
      ref = ref,
      within = within,
      method = method,
      test = if (method == "pairwise") test else NA,
      alpha = alpha,
      fdr_method = fdr_method
    ),
    data = data,
    call = call
  )

  validate_pepdiff_results(result_obj)
}


# =============================================================================
# GLM Method
# =============================================================================

#' Compare using GLM
#'
#' @keywords internal
compare_glm <- function(data, compare, ref, within, alpha, fdr_method) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required for GLM method", call. = FALSE)
  }

  # Run models for all peptides
  model_results <- run_models(data, compare, ref = ref, method = "glm")

  # Extract diagnostics
  diagnostics <- extract_diagnostics(model_results)

  # Build results tibble
  results <- build_results_tibble(
    model_results = model_results,
    data = data,
    compare = compare,
    method = "glm",
    alpha = alpha,
    fdr_method = fdr_method
  )

  list(results = results, diagnostics = diagnostics)
}


# =============================================================================
# ART Method
# =============================================================================

#' Compare using ART
#'
#' @keywords internal
compare_art <- function(data, compare, ref, within, alpha, fdr_method) {
  if (!requireNamespace("ARTool", quietly = TRUE)) {
    stop("Package 'ARTool' is required for ART method", call. = FALSE)
  }

  # Run models for all peptides
  model_results <- run_models(data, compare, ref = ref, method = "art")

  # Extract diagnostics
  diagnostics <- extract_diagnostics(model_results)

  # Build results tibble
  results <- build_results_tibble(
    model_results = model_results,
    data = data,
    compare = compare,
    method = "art",
    alpha = alpha,
    fdr_method = fdr_method
  )

  list(results = results, diagnostics = diagnostics)
}


# =============================================================================
# Pairwise Method
# =============================================================================

#' Compare using pairwise tests
#'
#' @keywords internal
compare_pairwise <- function(data, compare, ref, within, test, alpha, fdr_method) {
  peptides <- data$peptides

  # Get treatment level(s)
  all_levels <- unique(data$data[[compare]])
  treatment_levels <- setdiff(all_levels, ref)

  if (length(treatment_levels) > 1) {
    message("Note: pairwise method compares each treatment level to reference separately")
  }

  # Select the test function
  test_fn <- switch(test,
    "wilcoxon" = test_wilcoxon,
    "bootstrap_t" = test_bootstrap_t,
    "bayes_t" = test_bayes_t,
    "rankprod" = test_rankprod,
    stop("Unknown test: ", test, call. = FALSE)
  )

  all_results <- list()

  # For each treatment level (if multiple)
  for (trt_level in treatment_levels) {
    # Run test for each peptide
    test_results <- lapply(peptides, function(pep) {
      pep_data <- data$data[data$data$peptide == pep, ]

      # Get control and treatment values
      ctrl_data <- pep_data[pep_data[[compare]] == ref, ]
      trt_data <- pep_data[pep_data[[compare]] == trt_level, ]

      ctrl_values <- ctrl_data$value
      trt_values <- trt_data$value

      # Calculate fold change
      ctrl_mean <- mean(ctrl_values, na.rm = TRUE)
      trt_mean <- mean(trt_values, na.rm = TRUE)
      fc <- if (ctrl_mean > 0) trt_mean / ctrl_mean else NA_real_

      # Run test
      test_result <- test_fn(ctrl_values, trt_values)
      test_result$fold_change <- fc

      test_result
    })
    names(test_results) <- peptides

    # Build results for this comparison
    comparison_name <- paste(trt_level, "vs", ref)

    results_df <- tibble::tibble(
      peptide = peptides,
      gene_id = vapply(peptides, function(pep) {
        unique(data$data$gene_id[data$data$peptide == pep])[1]
      }, character(1)),
      comparison = comparison_name,
      fold_change = vapply(test_results, function(x) x$fold_change %||% NA_real_, numeric(1)),
      log2_fc = vapply(test_results, function(x) {
        fc <- x$fold_change %||% NA_real_
        if (is.na(fc) || fc <= 0) NA_real_ else log2(fc)
      }, numeric(1)),
      test = test,
      p_value = vapply(test_results, function(x) x$p_value %||% NA_real_, numeric(1))
    )

    all_results[[comparison_name]] <- results_df
  }

  # Combine all comparisons
  results <- dplyr::bind_rows(all_results)

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

  # Diagnostics for pairwise (simpler - just track which had valid results)
  diagnostics <- tibble::tibble(
    peptide = peptides,
    converged = !is.na(results$p_value[match(peptides, results$peptide)])
  )

  list(results = results, diagnostics = diagnostics)
}


# =============================================================================
# Stratified Analysis Helper
# =============================================================================

#' Compare within strata
#'
#' Runs comparisons separately within each level of stratifying factor(s).
#'
#' @param data pepdiff_data object
#' @param compare Factor to compare
#' @param ref Reference level
#' @param within Factor(s) to stratify by
#' @param method Analysis method
#' @param test Pairwise test (if applicable)
#' @param alpha Significance threshold
#' @param fdr_method FDR correction method
#'
#' @return Combined results across strata
#' @keywords internal
compare_within_strata <- function(data, compare, ref, within, method, test, alpha, fdr_method) {
  # Get all unique combinations of within factors
  within_data <- unique(data$data[, within, drop = FALSE])

  all_results <- list()

  for (i in seq_len(nrow(within_data))) {
    # Subset data for this stratum
    stratum <- within_data[i, ]
    stratum_name <- paste(sapply(names(stratum), function(n) paste0(n, "=", stratum[[n]])), collapse = ", ")

    # Filter data to this stratum
    subset_data <- data
    subset_idx <- rep(TRUE, nrow(data$data))
    for (w in within) {
      subset_idx <- subset_idx & (data$data[[w]] == stratum[[w]])
    }
    subset_data$data <- data$data[subset_idx, ]
    subset_data$peptides <- unique(subset_data$data$peptide)

    # Run analysis on subset
    if (method == "pairwise") {
      stratum_results <- compare_pairwise(subset_data, compare, ref, NULL, test, alpha, fdr_method)
    } else if (method == "glm") {
      stratum_results <- compare_glm(subset_data, compare, ref, NULL, alpha, fdr_method)
    } else {
      stratum_results <- compare_art(subset_data, compare, ref, NULL, alpha, fdr_method)
    }

    # Add stratum info to results
    stratum_results$results$stratum <- stratum_name
    for (w in within) {
      stratum_results$results[[w]] <- stratum[[w]]
    }

    all_results[[stratum_name]] <- stratum_results
  }

  # Combine results
  combined_results <- dplyr::bind_rows(lapply(all_results, function(x) x$results))
  combined_diagnostics <- dplyr::bind_rows(lapply(all_results, function(x) x$diagnostics))

  list(results = combined_results, diagnostics = combined_diagnostics)
}
