# Internal utility functions for pepdiff
# This file contains validation helpers and internal utilities

#' @importFrom rlang %||%
NULL

# =============================================================================
# Validation Functions
# =============================================================================

#' Validate factor columns exist in data
#'
#' Checks that all specified factor columns exist in the data frame
#' and are not empty.
#'
#' @param data A data frame
#' @param factors Character vector of factor column names
#' @return TRUE invisibly if valid, otherwise throws an error
#' @keywords internal
validate_factors <- function(data, factors) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  if (length(factors) == 0) {
    stop("At least one factor must be specified", call. = FALSE)
  }

  if (!is.character(factors)) {
    stop("'factors' must be a character vector", call. = FALSE)
  }

  missing <- setdiff(factors, names(data))
  if (length(missing) > 0) {
    stop(
      "Factor columns not found in data: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Check factors have at least 2 levels

for (f in factors) {
    unique_vals <- unique(data[[f]][!is.na(data[[f]])])
    if (length(unique_vals) < 2) {
      stop(
        "Factor '", f, "' must have at least 2 levels, found: ",
        length(unique_vals),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


#' Validate positive values for Gamma GLM
#'
#' Checks that all non-NA values are strictly positive (> 0).
#' Gamma GLM requires positive values.
#'
#' @param values Numeric vector of values to check
#' @param name Name of the variable for error messages
#' @return TRUE invisibly if valid, otherwise throws an error
#' @keywords internal
validate_positive <- function(values, name = "values") {
  if (!is.numeric(values)) {
    stop("'", name, "' must be numeric", call. = FALSE)
  }

  non_na_values <- values[!is.na(values)]

  if (any(non_na_values <= 0)) {
    n_zero <- sum(non_na_values == 0)
    n_negative <- sum(non_na_values < 0)
    msg <- paste0(
      "Gamma GLM requires strictly positive values. Found in '", name, "': "
    )
    parts <- c()
    if (n_zero > 0) parts <- c(parts, paste0(n_zero, " zeros"))
    if (n_negative > 0) parts <- c(parts, paste0(n_negative, " negative values"))
    stop(msg, paste(parts, collapse = " and "), call. = FALSE)
  }

  invisible(TRUE)
}


#' Check for zeros in values (used before GLM fitting)
#'
#' @param values Numeric vector
#' @return TRUE invisibly if no zeros/negatives, otherwise throws error
#' @keywords internal
check_zeros <- function(values) {
  validate_positive(values, "abundance values")
}


# =============================================================================
# Missingness Computation
# =============================================================================

#' Compute missingness statistics per peptide
#'
#' Calculates the NA rate and MNAR (Missing Not At Random) score for each peptide.
#' MNAR score is based on the relationship between missingness and mean abundance.
#'
#' @param data A data frame with 'peptide' and 'value' columns
#' @return A tibble with columns: peptide, na_rate, mnar_score, mean_abundance
#' @keywords internal
#' @importFrom rlang .data
compute_missingness <- function(data) {
  if (!all(c("peptide", "value") %in% names(data))) {
    stop("Data must have 'peptide' and 'value' columns", call. = FALSE)
  }

  stats <- data %>%
    dplyr::group_by(.data$peptide) %>%
    dplyr::summarize(
      n_total = dplyr::n(),
      n_missing = sum(is.na(.data$value)),
      mean_abundance = mean(.data$value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      na_rate = .data$n_missing / .data$n_total
    )

  # Compute MNAR score: correlation between missingness and low abundance
  # Higher score = more likely MNAR (missing values correlate with low abundance)
  # Using rank-based approach: peptides with more missing data and lower abundance
  # get higher MNAR scores
  if (nrow(stats) > 1 && any(stats$na_rate > 0)) {
    # Rank by abundance (lower = higher rank number)
    abundance_rank <- rank(-stats$mean_abundance, na.last = "keep")
    # MNAR score: weighted combination of na_rate and low abundance
    stats$mnar_score <- stats$na_rate * (abundance_rank / max(abundance_rank, na.rm = TRUE))
  } else {
    stats$mnar_score <- 0
  }

  stats %>%
    dplyr::select(
      "peptide",
      "na_rate",
      "mnar_score",
      "mean_abundance"
    )
}


# =============================================================================
# Design Summary
# =============================================================================

#' Compute design summary from data
#'
#' Creates a summary of factor combinations with replicate and peptide counts.
#'
#' @param data A data frame with factor columns, 'bio_rep', and 'peptide' columns
#' @param factors Character vector of factor column names
#' @return A tibble with factor columns plus n_reps and n_peptides
#' @keywords internal
#' @importFrom rlang .data
compute_design_summary <- function(data, factors) {
  if (!all(c("bio_rep", "peptide") %in% names(data))) {
    stop("Data must have 'bio_rep' and 'peptide' columns", call. = FALSE)
  }

  # Group by all factors
  group_vars <- c(factors)

  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarize(
      n_reps = dplyr::n_distinct(.data$bio_rep),
      n_peptides = dplyr::n_distinct(.data$peptide),
      n_observations = dplyr::n(),
      .groups = "drop"
    )
}


# =============================================================================
# Helper Functions
# =============================================================================

#' Check if values are useable (non-NA and non-NaN)
#'
#' @param x Vector of values
#' @return Logical vector indicating positions of useable values
#' @keywords internal
is_useable <- function(x) {
  !is.na(x) & !is.nan(x)
}


#' Calculate mean fold change for peptide
#'
#' For each peptide calculate the mean quantification over all experiments
#' then get the natural scale fold change.
#'
#' @param treatment Vector or matrix of treatment data
#' @param control Vector or matrix of control data
#' @return Vector of mean fold changes: mean(treatment) / mean(control)
#' @keywords internal
mean_fold_change <- function(treatment, control) {
  if (is.matrix(treatment)) {
    rowMeans(treatment, na.rm = TRUE) / rowMeans(control, na.rm = TRUE)
  } else {
    mean(treatment, na.rm = TRUE) / mean(control, na.rm = TRUE)
  }
}


#' Calculate log2 fold change
#'
#' @param fold_change Numeric vector of fold changes
#' @return Numeric vector of log2 fold changes
#' @keywords internal
log2_fc <- function(fold_change) {
  log2(fold_change)
}


#' Apply FDR correction within groups
#'
#' Applies Benjamini-Hochberg FDR correction within each comparison group.
#'
#' @param results A data frame with 'comparison' and 'p_value' columns
#' @param method FDR correction method (default "BH")
#' @return Data frame with added 'fdr' column
#' @keywords internal
#' @importFrom rlang .data
apply_fdr_by_comparison <- function(results, method = "BH") {
  results %>%
    dplyr::group_by(.data$comparison) %>%
    dplyr::mutate(
      fdr = stats::p.adjust(.data$p_value, method = method)
    ) %>%
    dplyr::ungroup()
}


#' Build formula from factor names
#'
#' Creates a formula string for GLM fitting.
#'
#' @param response Name of response variable
#' @param factors Character vector of factor names
#' @param interaction Logical, include interactions? (default TRUE)
#' @return A formula object
#' @keywords internal
build_formula <- function(response, factors, interaction = TRUE) {
  if (interaction && length(factors) > 1) {
    rhs <- paste(factors, collapse = " * ")
  } else {
    rhs <- paste(factors, collapse = " + ")
  }
  stats::as.formula(paste(response, "~", rhs))
}


# =============================================================================
# Legacy Compatibility Helpers (for internal use)
# =============================================================================

#' Convert wide format results table to long format
#'
#' Tidies up the wide results table from `compare()` to a long format.
#'
#' @param r Results dataframe typically from `compare()`
#' @return Dataframe in long format
#' @export
#' @importFrom rlang .data
long_results <- function(r) {
  r %>%
    dplyr::select(
      .data$gene_id,
      .data$peptide,
      dplyr::any_of(c("treatment_replicates", "control_replicates")),
      .data$fold_change,
      dplyr::ends_with("p_val")
    ) %>%
    tidyr::pivot_longer(dplyr::ends_with("p_val"), names_to = "test")
}
