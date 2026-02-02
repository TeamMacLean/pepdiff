# Legacy compatibility functions for pepdiff v1
# These functions wrap the old API with deprecation warnings

# =============================================================================
# Deprecated compare() default method
# =============================================================================

#' Default compare method for legacy data frames
#'
#' This method handles calls to compare() with data frames from the old
#' import_data() function. It issues a deprecation warning and delegates
#' to compare_legacy().
#'
#' @param data A data frame from import_data()
#' @param ... Arguments passed to compare_legacy()
#' @return Results from compare_legacy()
#' @export
#' @method compare data.frame
compare.data.frame <- function(data, ...) {
  .Deprecated("compare.pepdiff_data")
  message("Note: Use read_pepdiff() and compare.pepdiff_data() for new analyses")
  compare_legacy(data, ...)
}


# =============================================================================
# Legacy helper functions needed by compare_legacy
# =============================================================================

#' Convert data frame to matrix format (legacy)
#'
#' @param df Data frame from import_data()
#' @return List with row_info and data matrix
#' @keywords internal
#' @importFrom rlang .data
matrix_data <- function(df) {
  df <- combine_tech_reps(df)
  row_info <- dplyr::group_by(df, .data$gene_id, .data$peptide, .data$treatment, .data$seconds, .data$bio_rep) %>%
    dplyr::summarize(col_count = dplyr::n(), .groups = "drop")

  dm <- df %>%
    tidyr::pivot_wider(names_from = c(.data$treatment, .data$seconds, .data$bio_rep), values_from = .data$mean_tr_quant) %>%
    as.matrix()

  row_info <- dm[, c("gene_id", "peptide")]
  col_info <- colnames(dm[, 3:length(colnames(dm))])
  dm <- matrix(as.numeric(dm[, 3:ncol(dm)]), nrow = nrow(dm))
  colnames(dm) <- col_info
  return(list(row_info = row_info, data = dm))
}


#' Select columns for contrast (legacy)
#'
#' @param l List of data and row_info
#' @param treatment Treatment name
#' @param t_seconds Treatment seconds
#' @param control Control name
#' @param c_seconds Control seconds
#' @return List with treatment and control matrices
#' @keywords internal
select_columns_for_contrast <- function(l, treatment = NA,
                                        t_seconds = NA,
                                        control = NA,
                                        c_seconds = NA) {
  td <- paste(treatment, t_seconds, sep = "_")
  td <- paste0("^", td, "$")

  cd <- paste(control, c_seconds, sep = "_")
  cd <- paste0("^", cd, "$")
  cnames_no_biorep <- unlist(lapply(strsplit(colnames(l$data), "_"), function(x) paste0(x[1], "_", x[2])))

  t_ind <- which(stringr::str_detect(cnames_no_biorep, td))
  c_ind <- which(stringr::str_detect(cnames_no_biorep, cd))

  return(list(
    treatment = l$data[, t_ind],
    control = l$data[, c_ind]
  ))
}


#' Find minimal values for a peptide across replicates (legacy)
#'
#' @param d Matrix of peptide quantifications
#' @return Vector of minimum values
#' @keywords internal
min_peptide_values <- function(d) {
  apply(d$data, MARGIN = 1, min, na.rm = TRUE)
}


#' Replace missing values with replacements (legacy)
#'
#' @param x Vector with missing values
#' @param lowest_vals Vector of replacements
#' @return Vector with values replaced
#' @keywords internal
replace_vals <- function(x, lowest_vals) {
  to_replace <- which(is.na(x))
  x[to_replace] <- lowest_vals[to_replace]
  return(x)
}


#' Get p-values using normal percentile (legacy)
#'
#' @param treatment Treatment matrix
#' @param control Control matrix
#' @param iters Number of iterations
#' @return Data frame with p-values
#' @keywords internal
get_percentile_lowest_observed_value_iterative <- function(treatment, control, iters = 1000) {
  peptide_means <- apply(control, MARGIN = 1, mean, na.rm = TRUE)
  peptide_sds <- apply(control, MARGIN = 1, stats::sd, na.rm = TRUE)

  nobs <- dim(control)[1]
  m <- matrix(NA, nrow = nobs, ncol = iters)
  for (i in 1:iters) {
    m[, i] <- stats::rnorm(nobs, peptide_means, peptide_sds)
  }

  t_means <- apply(treatment, MARGIN = 1, mean, na.rm = TRUE)

  find_p <- function(x, perc) {
    stats::ecdf(x)(perc)
  }
  r <- rep(NA, nobs)
  for (i in 1:nobs) {
    r[i] <- find_p(m[i, ], t_means[i])
  }
  return(data.frame(norm_quantile_pval = r))
}
