# GLM Fit Diagnostics for pepdiff
# Provides visual diagnostics to assess GLM model fit across peptides

#' Plot GLM fit diagnostics
#'
#' Creates a multi-panel diagnostic plot to help assess whether GLM models
#' fit the data well. This is useful for deciding whether to use GLM or
#' switch to ART (Aligned Rank Transform).
#'
#' @param results A `pepdiff_results` object from `compare()` with `method = "glm"`
#' @param n_sample Number of peptides to show in sample residual plots (default 6)
#' @param deviance_threshold Optional threshold for flagging high-deviance peptides.
#'   If NULL (default), uses the 95th percentile of deviance values.
#' @param full_qq Deprecated. Residuals are now stored during `compare()` so
#'   accurate QQ plots are always available without refitting.
#'
#' @return Invisibly returns a list with:
#'   \item{plot}{The diagnostic plot (ggplot/cowplot grid)}
#'   \item{flagged}{Tibble of peptides with potential fit issues}
#'   \item{summary}{List with summary statistics (n_flagged, median_deviance, etc.)}
#'
#' @details
#' The function generates a 4-panel diagnostic plot:
#'
#' **Panel 1: Deviance Distribution** - Histogram showing the distribution of
#' residual deviance across all converged peptides. A long right tail suggests
#' some peptides fit poorly.
#'
#' **Panel 2: Deviance vs Fold Change** - Scatter plot of deviance against
#' absolute log2 fold change. If high-deviance points cluster at extreme
#' fold changes, this may indicate outlier-driven "significant" results.
#'
#' **Panel 3: Sample Residual Plots** - Residuals vs fitted values for a
#' sample of peptides (2 with highest deviance, 2 median, 2 lowest). Look
#' for random scatter around zero; patterns or funnels indicate poor fit.
#'
#' **Panel 4: Pooled QQ Plot** - Quantile-quantile plot of pooled residuals.
#' Points should fall on the diagonal line. S-curves indicate heavy tails
#' (consider ART), systematic deviation suggests wrong distributional assumption.
#'
#' @section Interpretation:
#'
#' **Use GLM when:**
#' - Deviance distribution looks reasonable (few flagged peptides)
#' - No systematic patterns in residual plots
#' - QQ plot is reasonably linear
#'
#' **Consider ART when:**
#' - Many peptides (>15%) have high deviance
#' - Residual plots show systematic curves or funnels
#' - QQ plot shows heavy tails (S-curve)
#'
#' @seealso [compare()] for running the analysis, `vignette("checking_fit")`
#'   for detailed guidance on interpreting diagnostics
#'
#' @examples
#' \dontrun{
#' # Run GLM analysis
#' results <- compare(dat, compare = "treatment", ref = "ctrl", method = "glm")
#'
#' # Check fit diagnostics
#' diag <- plot_fit_diagnostics(results)
#'
#' # View flagged peptides
#' diag$flagged
#'
#' # Get summary statistics
#' diag$summary
#' }
#'
#' @export
#' @importFrom rlang .data
plot_fit_diagnostics <- function(results, n_sample = 6, deviance_threshold = NULL, full_qq = FALSE) {
  # Validate input
  if (!inherits(results, "pepdiff_results")) {
    stop("Input must be a pepdiff_results object from compare()", call. = FALSE)
  }

  # Check method is GLM

  if (results$method != "glm") {
    stop(
      "ART is non-parametric; residual diagnostics don't apply.\n",
      "See vignette('checking_fit') for guidance on when to use ART.",
      call. = FALSE
    )
  }

  # Extract diagnostics with deviance values
  diag_data <- results$diagnostics
  if (is.null(diag_data) || !"deviance" %in% names(diag_data)) {
    stop("Diagnostics data not available. Re-run compare() with method = 'glm'.", call. = FALSE)
  }

  # Filter to converged peptides with valid deviance
  diag_data <- diag_data[diag_data$converged & !is.na(diag_data$deviance), ]

  if (nrow(diag_data) == 0) {
    stop("No converged peptides with valid deviance values.", call. = FALSE)
  }

  # Join with results to get fold changes and p-values
  plot_data <- dplyr::left_join(
    diag_data,
    results$results[, c("peptide", "gene_id", "fold_change", "log2_fc", "p_value", "significant")],
    by = "peptide"
  )

  # Handle duplicates from multiple comparisons (take first)
  plot_data <- plot_data[!duplicated(plot_data$peptide), ]

  # Determine threshold
  if (is.null(deviance_threshold)) {
    deviance_threshold <- stats::quantile(plot_data$deviance, 0.95, na.rm = TRUE)
  }

  # Calculate summary statistics
  n_analyzed <- nrow(plot_data)
  median_deviance <- stats::median(plot_data$deviance, na.rm = TRUE)
  n_flagged <- sum(plot_data$deviance > deviance_threshold, na.rm = TRUE)
  pct_flagged <- round(100 * n_flagged / n_analyzed, 1)

  # Create flagged peptides tibble
  flagged <- plot_data[plot_data$deviance > deviance_threshold, ]
  flagged <- tibble::tibble(
    peptide = flagged$peptide,
    gene_id = flagged$gene_id,
    deviance = flagged$deviance,
    fold_change = flagged$fold_change,
    p_value = flagged$p_value,
    significant = flagged$significant,
    flag_reason = "high_deviance"
  )

  # Build the 4-panel plot
  p1 <- build_deviance_histogram(plot_data, deviance_threshold, n_flagged)
  p2 <- build_deviance_vs_fc(plot_data, deviance_threshold)
  p3 <- build_sample_residual_plots(results, plot_data, n_sample)
  p4 <- build_pooled_qq(results, plot_data, full_qq)

  # Combine into grid
  combined_plot <- cowplot::plot_grid(
    p1, p2, p3, p4,
    ncol = 2,
    labels = c("A", "B", "C", "D")
  )

  # Build summary list
  summary_stats <- list(
    n_analyzed = n_analyzed,
    n_converged = sum(results$diagnostics$converged, na.rm = TRUE),
    n_failed = sum(!results$diagnostics$converged, na.rm = TRUE),
    median_deviance = median_deviance,
    threshold = deviance_threshold,
    n_flagged = n_flagged,
    pct_flagged = pct_flagged
  )

  # Print summary to console
  print_diagnostic_summary(summary_stats)

  # Return invisibly
  invisible(list(
    plot = combined_plot,
    flagged = flagged,
    summary = summary_stats
  ))
}


# =============================================================================
# Panel 1: Deviance Distribution Histogram
# =============================================================================

#' @keywords internal
build_deviance_histogram <- function(plot_data, threshold, n_flagged) {
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$deviance)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::labs(
      title = "Deviance Distribution",
      subtitle = paste(n_flagged, "peptides above threshold (red line)"),
      x = "Residual Deviance",
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40")
    )
}


# =============================================================================
# Panel 2: Deviance vs Fold Change
# =============================================================================

#' @keywords internal
build_deviance_vs_fc <- function(plot_data, threshold) {
  # Filter out NA fold changes
  fc_data <- plot_data[!is.na(plot_data$log2_fc), ]

  if (nrow(fc_data) == 0) {
    # Return empty plot with message
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No fold change data available") +
        ggplot2::theme_void()
    )
  }

  ggplot2::ggplot(fc_data, ggplot2::aes(x = abs(.data$log2_fc), y = .data$deviance)) +
    ggplot2::geom_point(alpha = 0.5, color = "steelblue") +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::labs(
      title = "Deviance vs Effect Size",
      subtitle = "High deviance at extreme FC may indicate outlier-driven results",
      x = "|log2 Fold Change|",
      y = "Deviance"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40")
    )
}


# =============================================================================
# Panel 3: Sample Residual Plots
# =============================================================================

#' @keywords internal
build_sample_residual_plots <- function(results, plot_data, n_sample) {
  # Select peptides: high deviance, median, low deviance
  n_each <- floor(n_sample / 3)
  if (n_each < 1) n_each <- 1

  ordered_peps <- plot_data$peptide[order(plot_data$deviance, decreasing = TRUE)]

  # High deviance (top)
  high_peps <- utils::head(ordered_peps, n_each)

  # Low deviance (bottom)
  low_peps <- utils::tail(ordered_peps, n_each)

  # Median deviance (middle)
  n_total <- length(ordered_peps)
  mid_start <- max(1, floor(n_total / 2) - floor(n_each / 2))
  mid_peps <- ordered_peps[mid_start:min(mid_start + n_each - 1, n_total)]

  selected_peps <- unique(c(high_peps, mid_peps, low_peps))

  # Extract stored residuals and fitted values
  diag <- results$diagnostics
  resid_data <- extract_stored_residuals(diag, selected_peps)

  if (nrow(resid_data) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Could not extract residuals") +
        ggplot2::theme_void()
    )
  }

  # Order peptides by deviance category for display
  pep_deviance <- plot_data$deviance[match(resid_data$peptide, plot_data$peptide)]
  resid_data$deviance_rank <- factor(
    dplyr::case_when(
      resid_data$peptide %in% high_peps ~ "1_High",
      resid_data$peptide %in% mid_peps ~ "2_Median",
      TRUE ~ "3_Low"
    )
  )

  # Shorten peptide names for display
  resid_data$pep_label <- paste0(
    gsub("PEP_", "", resid_data$peptide),
    " (",
    resid_data$deviance_rank,
    ")"
  )

  ggplot2::ggplot(resid_data, ggplot2::aes(x = .data$fitted, y = .data$residual)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8) +
    ggplot2::facet_wrap(~ .data$peptide, scales = "free", ncol = 3) +
    ggplot2::labs(
      title = "Sample Residual Plots",
      subtitle = "High/median/low deviance peptides | Flat red line = good fit",
      x = "Fitted Values",
      y = "Deviance Residuals"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40"),
      strip.text = ggplot2::element_text(size = 8)
    )
}


#' Extract stored residuals and fitted values for selected peptides
#' @keywords internal
extract_stored_residuals <- function(diagnostics, peptides) {

  resid_list <- lapply(peptides, function(pep) {
    idx <- which(diagnostics$peptide == pep)
    if (length(idx) == 0) return(NULL)

    # Use single bracket indexing for tibble list columns, then extract element
    resid <- diagnostics$residuals[idx][[1]]
    fitted <- diagnostics$fitted[idx][[1]]

    if (is.null(resid) || is.null(fitted)) return(NULL)

    data.frame(
      peptide = pep,
      fitted = fitted,
      residual = resid
    )
  })

  dplyr::bind_rows(resid_list)
}


# =============================================================================
# Panel 4: Pooled QQ Plot
# =============================================================================

#' @keywords internal
build_pooled_qq <- function(results, plot_data, full_qq) {
  # Use stored standardized residuals from diagnostics (for QQ plot)
  diag <- results$diagnostics

  if ("std_residuals" %in% names(diag)) {
    # Use stored standardized residuals (preferred for QQ plot)
    pooled_resid <- unlist(diag$std_residuals[diag$converged])
  } else if ("residuals" %in% names(diag)) {
    # Fallback to raw residuals (older results objects)
    pooled_resid <- unlist(diag$residuals[diag$converged])
  } else {
    # Last resort: refit models
    pooled_resid <- get_all_residuals(results)
  }

  if (length(pooled_resid) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Could not compute residuals") +
        ggplot2::theme_void()
    )
  }

  qq_data <- data.frame(resid = pooled_resid)

  ggplot2::ggplot(qq_data, ggplot2::aes(sample = .data$resid)) +
    ggplot2::stat_qq(alpha = 0.5, color = "steelblue") +
    ggplot2::stat_qq_line(color = "red", linewidth = 1) +
    ggplot2::labs(
      title = "Pooled Residuals QQ Plot",
      subtitle = "Standardized deviance residuals from all peptides",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9, color = "gray40")
    )
}


#' Get standardized residuals from all peptide models
#' @keywords internal
get_all_residuals <- function(results) {
  data <- results$data
  factors <- data$factors
  peptides <- data$peptides

  formula <- build_formula("value", factors, interaction = length(factors) > 1)

  all_resid <- lapply(peptides, function(pep) {
    pep_data <- data$data[data$data$peptide == pep, ]
    pep_data <- pep_data[!is.na(pep_data$value), ]

    if (nrow(pep_data) < 3) return(NULL)
    if (any(pep_data$value <= 0)) return(NULL)

    tryCatch({
      model <- stats::glm(
        formula,
        data = pep_data,
        family = stats::Gamma(link = "log")
      )

      if (!model$converged) return(NULL)

      # Use standardized deviance residuals
      stats::rstandard(model, type = "deviance")
    }, error = function(e) NULL)
  })

  unlist(all_resid)
}


# =============================================================================
# Console Output
# =============================================================================

#' Print diagnostic summary to console
#' @keywords internal
print_diagnostic_summary <- function(summary_stats) {
  cat("\n")
  cat("GLM Fit Diagnostics\n")
  cat("-------------------\n")
  cat(sprintf("Peptides analyzed: %d", summary_stats$n_analyzed))
  if (summary_stats$n_failed > 0) {
    cat(sprintf(" (%d failed to converge)", summary_stats$n_failed))
  } else {
    cat(" (all converged)")
  }
  cat("\n")
  cat(sprintf("Median deviance: %.2f\n", summary_stats$median_deviance))
  cat(sprintf("Threshold: %.2f (95th percentile)\n", summary_stats$threshold))
  cat(sprintf("Flagged peptides: %d (%.1f%%) above threshold\n",
              summary_stats$n_flagged, summary_stats$pct_flagged))
  cat("\n")

  # Interpretation guidance
  if (summary_stats$pct_flagged < 5) {
    cat("Interpretation: Deviance distribution looks good.\n")
  } else if (summary_stats$pct_flagged < 15) {
    cat("Interpretation: Some peptides show elevated deviance. Check residual plots.\n")
  } else {
    cat("Interpretation: Many peptides show poor fit. Consider using method = 'art'.\n")
  }

  cat("See vignette('checking_fit') for detailed guidance.\n")
  cat("\n")
}
