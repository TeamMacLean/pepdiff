# GLM and ART model fitting for differential abundance analysis
# These functions handle Gamma GLM and Aligned Rank Transform methods

# =============================================================================
# GLM Fitting (Gamma with log link)
# =============================================================================

#' Fit a Gamma GLM for a single peptide
#'
#' Fits a generalized linear model with Gamma family and log link to model
#' abundance values as a function of experimental factors.
#'
#' @param data Data frame containing the peptide data
#' @param formula A formula specifying the model
#' @param peptide_id Character, the peptide being analyzed (for diagnostics)
#'
#' @return A list with components:
#'   \item{converged}{Logical, whether the model converged}
#'   \item{model}{The fitted glm object (NULL if failed)}
#'   \item{coefficients}{Model coefficients (NULL if failed)}
#'   \item{deviance}{Residual deviance (NA if failed)}
#'   \item{peptide}{The peptide ID}
#'
#' @keywords internal
fit_glm <- function(data, formula, peptide_id = NULL) {
  # Check for zeros (Gamma GLM requires positive values)
  values <- data[[all.vars(formula)[1]]]
  values <- values[!is.na(values)]

  if (any(values <= 0)) {
    return(list(
      converged = FALSE,
      model = NULL,
      coefficients = NULL,
      deviance = NA_real_,
      peptide = peptide_id,
      error = "Gamma GLM requires strictly positive values (zeros or negatives found)"
    ))
  }

  # Attempt to fit the model
  result <- tryCatch(
    {
      model <- stats::glm(
        formula,
        data = data,
        family = stats::Gamma(link = "log"),
        control = stats::glm.control(maxit = 100)
      )

      list(
        converged = model$converged,
        model = model,
        coefficients = stats::coef(model),
        deviance = model$deviance,
        peptide = peptide_id,
        error = if (!model$converged) "Model did not converge" else NULL
      )
    },
    warning = function(w) {
      # Try to continue despite warnings
      model <- suppressWarnings(stats::glm(
        formula,
        data = data,
        family = stats::Gamma(link = "log"),
        control = stats::glm.control(maxit = 100)
      ))

      list(
        converged = model$converged,
        model = model,
        coefficients = stats::coef(model),
        deviance = model$deviance,
        peptide = peptide_id,
        error = paste("Warning:", conditionMessage(w))
      )
    },
    error = function(e) {
      list(
        converged = FALSE,
        model = NULL,
        coefficients = NULL,
        deviance = NA_real_,
        peptide = peptide_id,
        error = conditionMessage(e)
      )
    }
  )

  result
}


#' Extract contrasts from a fitted GLM using emmeans
#'
#' Uses the emmeans package to compute estimated marginal means and
#' contrasts between factor levels.
#'
#' @param model A fitted glm object
#' @param specs Character vector of factor names to compare, or emmeans formula
#' @param ref Reference level for pairwise comparisons
#' @param adjust P-value adjustment method (default "none" - FDR done later)
#'
#' @return A tibble with columns:
#'   \item{contrast}{Description of the contrast}
#'   \item{estimate}{Log fold change estimate}
#'   \item{se}{Standard error}
#'   \item{z_ratio}{Z statistic}
#'   \item{p_value}{P-value}
#'   \item{fold_change}{Back-transformed fold change}
#'
#' @keywords internal
extract_contrasts_glm <- function(model, specs, ref = NULL, adjust = "none") {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("Package 'emmeans' is required for GLM contrasts", call. = FALSE)
  }

  # Get estimated marginal means
  emm <- emmeans::emmeans(model, specs = specs)

  # Compute pairwise contrasts
  if (!is.null(ref)) {
    # Treatment vs reference contrasts
    contrasts <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ref, adjust = adjust)
  } else {
    # All pairwise contrasts
    contrasts <- emmeans::contrast(emm, method = "pairwise", adjust = adjust)
  }

  # Convert to data frame
  contrast_df <- as.data.frame(contrasts)

  # Standardize column names
  result <- tibble::tibble(
    contrast = contrast_df$contrast,
    estimate = contrast_df$estimate,
    se = contrast_df$SE,
    z_ratio = contrast_df$z.ratio,
    p_value = contrast_df$p.value,
    # Back-transform from log scale to fold change
    fold_change = exp(contrast_df$estimate)
  )

  result
}


#' Fit GLM and extract contrasts for one peptide
#'
#' Convenience function that fits a GLM and extracts contrasts in one step.
#'
#' @param data Data for one peptide
#' @param response Name of the response variable
#' @param factors Character vector of factor names
#' @param compare Which factor to compare
#' @param ref Reference level
#' @param peptide_id Peptide identifier
#'
#' @return A list with model fit results and contrasts
#' @keywords internal
fit_and_extract_glm <- function(data, response, factors, compare, ref, peptide_id) {
  # Build formula
  formula <- build_formula(response, factors, interaction = length(factors) > 1)

  # Fit model
  fit_result <- fit_glm(data, formula, peptide_id)

  if (!fit_result$converged || is.null(fit_result$model)) {
    return(list(
      peptide = peptide_id,
      converged = FALSE,
      contrasts = NULL,
      diagnostics = fit_result
    ))
  }

  # Extract contrasts
  contrasts <- tryCatch(
    {
      extract_contrasts_glm(fit_result$model, specs = compare, ref = ref)
    },
    error = function(e) {
      NULL
    }
  )

  list(
    peptide = peptide_id,
    converged = TRUE,
    contrasts = contrasts,
    diagnostics = fit_result
  )
}


# =============================================================================
# ART Fitting (Aligned Rank Transform)
# =============================================================================

#' Fit an Aligned Rank Transform model for a single peptide
#'
#' Fits a non-parametric ART model which is suitable for factorial designs
#' when parametric assumptions are not met.
#'
#' @param data Data frame containing the peptide data
#' @param formula A formula specifying the model
#' @param peptide_id Character, the peptide being analyzed
#'
#' @return A list with components:
#'   \item{converged}{Logical, whether fitting succeeded}
#'   \item{model}{The fitted art object (NULL if failed)}
#'   \item{peptide}{The peptide ID}
#'   \item{error}{Error message if failed}
#'
#' @keywords internal
fit_art <- function(data, formula, peptide_id = NULL) {
  if (!requireNamespace("ARTool", quietly = TRUE)) {
    return(list(
      converged = FALSE,
      model = NULL,
      peptide = peptide_id,
      error = "Package 'ARTool' is required for ART analysis"
    ))
  }

  result <- tryCatch(
    {
      model <- ARTool::art(formula, data = data)

      list(
        converged = TRUE,
        model = model,
        peptide = peptide_id,
        error = NULL
      )
    },
    error = function(e) {
      list(
        converged = FALSE,
        model = NULL,
        peptide = peptide_id,
        error = conditionMessage(e)
      )
    }
  )

  result
}


#' Extract contrasts from a fitted ART model
#'
#' Uses art.con() from ARTool to compute contrasts.
#'
#' @param model A fitted art object
#' @param specs Factor name or formula for contrasts
#' @param adjust P-value adjustment method
#'
#' @return A tibble with contrast results
#' @keywords internal
extract_contrasts_art <- function(model, specs, adjust = "none") {
  if (!requireNamespace("ARTool", quietly = TRUE)) {
    stop("Package 'ARTool' is required for ART contrasts", call. = FALSE)
  }

  # Get contrasts using art.con
  contrasts <- tryCatch(
    {
      ARTool::art.con(model, specs, adjust = adjust)
    },
    error = function(e) {
      return(NULL)
    }
  )

  if (is.null(contrasts)) {
    return(NULL)
  }

  # Convert to data frame
  contrast_df <- as.data.frame(contrasts)

  # Standardize column names (ARTool output varies)
  result <- tibble::tibble(
    contrast = if ("contrast" %in% names(contrast_df)) contrast_df$contrast else rownames(contrast_df),
    estimate = if ("estimate" %in% names(contrast_df)) contrast_df$estimate else NA_real_,
    se = if ("SE" %in% names(contrast_df)) contrast_df$SE else NA_real_,
    t_ratio = if ("t.ratio" %in% names(contrast_df)) contrast_df$t.ratio else NA_real_,
    p_value = if ("p.value" %in% names(contrast_df)) contrast_df$p.value else contrast_df[[ncol(contrast_df)]]
  )

  result
}


#' Fit ART and extract contrasts for one peptide
#'
#' @param data Data for one peptide
#' @param response Name of the response variable
#' @param factors Character vector of factor names
#' @param compare Which factor to compare
#' @param peptide_id Peptide identifier
#'
#' @return A list with model fit results and contrasts
#' @keywords internal
fit_and_extract_art <- function(data, response, factors, compare, peptide_id) {
  # Build formula (ART uses Error term for subject/replicate)
  if (length(factors) > 1) {
    formula <- stats::as.formula(paste(
      response, "~",
      paste(factors, collapse = " * ")
    ))
  } else {
    formula <- stats::as.formula(paste(response, "~", factors))
  }

  # Fit model
  fit_result <- fit_art(data, formula, peptide_id)

  if (!fit_result$converged || is.null(fit_result$model)) {
    return(list(
      peptide = peptide_id,
      converged = FALSE,
      contrasts = NULL,
      diagnostics = fit_result
    ))
  }

  # Extract contrasts
  contrasts <- tryCatch(
    {
      extract_contrasts_art(fit_result$model, specs = compare)
    },
    error = function(e) {
      NULL
    }
  )

  list(
    peptide = peptide_id,
    converged = TRUE,
    contrasts = contrasts,
    diagnostics = fit_result
  )
}


# =============================================================================
# Model Selection and Running
# =============================================================================

#' Run statistical model for all peptides
#'
#' Applies GLM or ART model to each peptide in the dataset.
#'
#' @param data A pepdiff_data object
#' @param compare Factor to compare
#' @param ref Reference level
#' @param method One of "glm" or "art"
#'
#' @return A list with results for each peptide
#' @keywords internal
run_models <- function(data, compare, ref = NULL, method = "glm") {
  if (!inherits(data, "pepdiff_data")) {
    stop("Data must be a pepdiff_data object", call. = FALSE)
  }

  peptides <- data$peptides
  factors <- data$factors

  results <- lapply(peptides, function(pep) {
    pep_data <- data$data[data$data$peptide == pep, ]

    # Remove rows with NA values
    pep_data <- pep_data[!is.na(pep_data$value), ]

    if (nrow(pep_data) < 3) {
      return(list(
        peptide = pep,
        converged = FALSE,
        contrasts = NULL,
        diagnostics = list(error = "Insufficient data points")
      ))
    }

    if (method == "glm") {
      fit_and_extract_glm(pep_data, "value", factors, compare, ref, pep)
    } else if (method == "art") {
      fit_and_extract_art(pep_data, "value", factors, compare, pep)
    } else {
      stop("Unknown method: ", method, call. = FALSE)
    }
  })

  names(results) <- peptides
  results
}


# =============================================================================
# Diagnostic Helpers
# =============================================================================

#' Extract model diagnostics
#'
#' Creates a summary tibble of model convergence and quality metrics.
#'
#' @param model_results List of model results from run_models()
#'
#' @return A tibble with diagnostic information
#' @keywords internal
extract_diagnostics <- function(model_results) {
  tibble::tibble(
    peptide = names(model_results),
    converged = vapply(model_results, function(x) x$converged, logical(1)),
    error = vapply(model_results, function(x) {
      err <- x$diagnostics$error
      if (is.null(err)) NA_character_ else err
    }, character(1)),
    deviance = vapply(model_results, function(x) {
      dev <- x$diagnostics$deviance
      if (is.null(dev)) NA_real_ else dev
    }, numeric(1))
  )
}


#' Check if a model has enough variation to fit
#'
#' @param values Numeric vector of values
#' @return Logical
#' @keywords internal
has_sufficient_variation <- function(values) {
  values <- values[!is.na(values)]
  if (length(values) < 3) return(FALSE)
  stats::var(values) > .Machine$double.eps
}
