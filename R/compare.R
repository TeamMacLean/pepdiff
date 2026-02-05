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
#' @param alpha Significance threshold (default 0.05). Used for p-value based tests.
#' @param fdr_method FDR correction method (default "BH"). Not applied for bayes_t.
#' @param bf_threshold Bayes factor threshold for significance (default 3).
#'   Only used when test = "bayes_t". BF > threshold marks peptide as significant.
#' @param ... Additional arguments passed to methods
#'
#' @return A pepdiff_results object containing:
#'   \item{results}{Tibble with peptide, gene_id, comparison, fold_change, log2_fc,
#'     p_value, fdr, significant. For bayes_t: p_value/fdr are NA, includes bf and evidence columns.}
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
#'
#' # Bayes factor test (uses bf_threshold instead of alpha/FDR)
#' results <- compare(data, compare = "treatment", ref = "ctrl",
#'                    method = "pairwise", test = "bayes_t", bf_threshold = 10)
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
                                  bf_threshold = 3,
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
  # If 'within' is specified, run stratified analysis for all methods
  if (!is.null(within)) {
    results <- compare_within_strata(data, compare, ref, within, method, test, alpha, fdr_method, bf_threshold)
  } else if (method == "pairwise") {
    results <- compare_pairwise(data, compare, ref, within, test, alpha, fdr_method, bf_threshold)
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
      fdr_method = fdr_method,
      bf_threshold = if (method == "pairwise" && test == "bayes_t") bf_threshold else NA
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
compare_pairwise <- function(data, compare, ref, within, test, alpha, fdr_method, bf_threshold = 3) {
  peptides <- data$peptides
  is_bayes <- test == "bayes_t"

  # Handle rankprod specially - needs all peptides together
  if (test == "rankprod") {
    return(compare_pairwise_rankprod(data, compare, ref, alpha, fdr_method))
  }

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

    if (is_bayes) {
      # Bayes factor results: bf and evidence columns, NA for p_value/fdr
      bf_values <- vapply(test_results, function(x) x$bf %||% NA_real_, numeric(1))

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
        p_value = NA_real_,
        fdr = NA_real_,
        bf = bf_values,
        evidence = classify_bf_evidence(bf_values),
        significant = bf_values > bf_threshold
      )
    } else {
      # P-value based tests
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
    }

    all_results[[comparison_name]] <- results_df
  }

  # Combine all comparisons
  results <- dplyr::bind_rows(all_results)

  if (!is_bayes) {
    # Apply FDR correction within each comparison (only for p-value tests)
    results <- results %>%
      dplyr::group_by(.data$comparison) %>%
      dplyr::mutate(
        fdr = stats::p.adjust(.data$p_value, method = fdr_method)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        significant = .data$fdr < alpha
      )
  }

  # Diagnostics for pairwise (simpler - just track which had valid results)
  if (is_bayes) {
    diagnostics <- tibble::tibble(
      peptide = peptides,
      converged = !is.na(results$bf[match(peptides, results$peptide)])
    )
  } else {
    diagnostics <- tibble::tibble(
      peptide = peptides,
      converged = !is.na(results$p_value[match(peptides, results$peptide)])
    )
  }

  list(results = results, diagnostics = diagnostics)
}


#' Compare using Rank Products (full matrix approach)
#'
#' Rank Products requires ranking across ALL peptides simultaneously,
#' not within single peptides. This function handles the matrix-based
#' approach using the RankProd package.
#'
#' @param data pepdiff_data object
#' @param compare Factor to compare
#' @param ref Reference level
#' @param alpha Significance threshold
#' @param fdr_method FDR correction method
#' @return List with results and diagnostics tibbles
#' @keywords internal
compare_pairwise_rankprod <- function(data, compare, ref, alpha, fdr_method) {
  peptides <- data$peptides

  # Check for RankProd package
  has_rankprod <- requireNamespace("RankProd", quietly = TRUE)
  if (!has_rankprod) {
    warning("RankProd package required for rankprod test. Install with BiocManager::install('RankProd'). Returning NA p-values.",
            call. = FALSE)
  }

  # Get treatment level(s)
  all_levels <- unique(data$data[[compare]])
  treatment_levels <- setdiff(all_levels, ref)

  if (length(treatment_levels) > 1) {
    message("Note: pairwise method compares each treatment level to reference separately")
  }

  all_results <- list()

  for (trt_level in treatment_levels) {
    comparison_name <- paste(trt_level, "vs", ref)

    if (!has_rankprod) {
      # Return NA results when package not available
      results_df <- tibble::tibble(
        peptide = peptides,
        gene_id = vapply(peptides, function(pep) {
          unique(data$data$gene_id[data$data$peptide == pep])[1]
        }, character(1)),
        comparison = comparison_name,
        fold_change = NA_real_,
        log2_fc = NA_real_,
        test = "rankprod",
        p_value = NA_real_,
        fdr = NA_real_,
        significant = NA
      )
      all_results[[comparison_name]] <- results_df
      next
    }

    # Build matrices (peptides × replicates)
    ctrl_mat <- build_replicate_matrix(data, compare, ref)
    trt_mat <- build_replicate_matrix(data, compare, trt_level)

    # Ensure consistent peptide order
    ctrl_mat <- ctrl_mat[peptides, , drop = FALSE]
    trt_mat <- trt_mat[peptides, , drop = FALSE]

    # RankProd expects: combined matrix with class labels
    # cl: 1 = treatment, 0 = control
    cl <- c(rep(1, ncol(trt_mat)), rep(0, ncol(ctrl_mat)))
    combined <- cbind(trt_mat, ctrl_mat)

    # Call RankProd (suppress its plot output)
    rp_result <- RankProd::RankProducts(
      combined, cl,
      logged = FALSE,
      plot = FALSE,
      na.rm = TRUE
    )

    # Extract p-values
    # Column 1: p-value for up-regulation (class 1 > class 0, i.e., treatment > control)
    # Column 2: p-value for down-regulation (class 1 < class 0)
    p_up <- rp_result$pval[, 1]
    p_down <- rp_result$pval[, 2]

    # Two-sided p-value: minimum of up/down
    p_value <- pmin(p_up, p_down)

    # Calculate fold changes
    ctrl_means <- rowMeans(ctrl_mat, na.rm = TRUE)
    trt_means <- rowMeans(trt_mat, na.rm = TRUE)
    fold_change <- ifelse(ctrl_means > 0, trt_means / ctrl_means, NA_real_)
    log2_fc <- ifelse(!is.na(fold_change) & fold_change > 0, log2(fold_change), NA_real_)

    results_df <- tibble::tibble(
      peptide = peptides,
      gene_id = vapply(peptides, function(pep) {
        unique(data$data$gene_id[data$data$peptide == pep])[1]
      }, character(1)),
      comparison = comparison_name,
      fold_change = fold_change,
      log2_fc = log2_fc,
      test = "rankprod",
      p_value = p_value
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

  # Diagnostics
  diagnostics <- tibble::tibble(
    peptide = peptides,
    converged = !is.na(results$p_value[match(peptides, results$peptide)])
  )

  list(results = results, diagnostics = diagnostics)
}


#' Build replicate matrix from pepdiff_data
#'
#' Pivots data to matrix format (peptides × replicates) for a given factor level.
#'
#' @param data pepdiff_data object
#' @param compare Factor column name
#' @param level Level of the factor to extract
#' @return Matrix with peptides as rows and replicates as columns
#' @keywords internal
build_replicate_matrix <- function(data, compare, level) {
  # Subset to the specified factor level
  subset_data <- data$data[data$data[[compare]] == level, ]

  # Pivot to wide format
  wide <- tidyr::pivot_wider(
    subset_data,
    id_cols = "peptide",
    names_from = "bio_rep",
    values_from = "value"
  )

  # Convert to matrix with peptide rownames
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$peptide
  mat
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
compare_within_strata <- function(data, compare, ref, within, method, test, alpha, fdr_method, bf_threshold = 3) {
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
    # Exclude 'within' factors from model (they're constant within stratum)
    subset_data$factors <- setdiff(data$factors, within)

    # Run analysis on subset
    if (method == "pairwise") {
      stratum_results <- compare_pairwise(subset_data, compare, ref, NULL, test, alpha, fdr_method, bf_threshold)
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
