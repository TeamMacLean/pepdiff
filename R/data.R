# pepdiff_data S3 class and data import functions
# Handles data import, preprocessing, and the core data container

# =============================================================================
# pepdiff_data Class Constructor and Validator
# =============================================================================

#' Create a new pepdiff_data object
#'
#' Low-level constructor for pepdiff_data objects. Use [read_pepdiff()]
#' for user-facing data import.
#'
#' @param data A tibble with peptide, gene_id, factor columns, bio_rep, value
#' @param factors Character vector of factor column names
#' @param design Tibble of unique factor combinations with counts
#' @param missingness Tibble of peptide missingness statistics
#' @param peptides Character vector of unique peptide IDs
#' @param call The original function call
#' @return A pepdiff_data object
#' @keywords internal
new_pepdiff_data <- function(data, factors, design, missingness, peptides, call) {
  structure(
    list(
      data = data,
      factors = factors,
      design = design,
      missingness = missingness,
      peptides = peptides,
      call = call
    ),
    class = "pepdiff_data"
  )
}


#' Validate a pepdiff_data object
#'
#' Checks that all required components are present and correctly structured.
#'
#' @param x A pepdiff_data object to validate
#' @return The validated object (invisibly), or throws an error
#' @keywords internal
validate_pepdiff_data <- function(x) {
  if (!inherits(x, "pepdiff_data")) {
    stop("Object must be of class 'pepdiff_data'", call. = FALSE)
  }

  # Check required components exist
  required <- c("data", "factors", "design", "missingness", "peptides", "call")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(
      "pepdiff_data missing required components: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Check data is a tibble/data.frame
  if (!is.data.frame(x$data)) {
    stop("'data' component must be a data frame", call. = FALSE)
  }

  # Check required columns in data
  required_cols <- c("peptide", "gene_id", "bio_rep", "value")
  missing_cols <- setdiff(required_cols, names(x$data))
  if (length(missing_cols) > 0) {
    stop(
      "Data missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Check factor columns exist in data
  missing_factors <- setdiff(x$factors, names(x$data))
  if (length(missing_factors) > 0) {
    stop(
      "Factor columns not found in data: ",
      paste(missing_factors, collapse = ", "),
      call. = FALSE
    )
  }

  # Check factors is character vector
  if (!is.character(x$factors) || length(x$factors) == 0) {
    stop("'factors' must be a non-empty character vector", call. = FALSE)
  }

  # Check peptides is character vector
  if (!is.character(x$peptides)) {
    stop("'peptides' must be a character vector", call. = FALSE)
  }

  invisible(x)
}


# =============================================================================
# Main Import Function
# =============================================================================

#' Read proteomics data into a pepdiff_data object
#'
#' Imports a CSV file containing PRM proteomics data and creates a
#' pepdiff_data object suitable for analysis with [compare()].
#'
#' @param file Path to CSV file
#' @param id Column name containing peptide identifiers
#' @param gene Column name containing gene identifiers
#' @param value Column name containing abundance values
#' @param factors Character vector of column names to use as experimental factors
#' @param replicate Column name containing biological replicate identifiers
#' @param tech_rep Optional column name containing technical replicate identifiers.
#'   If provided, data will NOT be automatically combined - use [combine_tech_reps()]
#'   explicitly after import.
#'
#' @return A pepdiff_data object with components:
#'   \item{data}{Tibble with columns: peptide, gene_id, [factors], bio_rep, value}
#'   \item{factors}{Character vector of factor names}
#'   \item{design}{Tibble of factor combinations with n_reps, n_peptides}
#'   \item{missingness}{Tibble of peptide missingness statistics}
#'   \item{peptides}{Character vector of unique peptide IDs}
#'   \item{call}{The original function call}
#'
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Simple import with one factor
#' dat <- read_pepdiff(
#'   "data.csv",
#'   id = "peptide_sequence",
#'   gene = "gene_name",
#'   value = "intensity",
#'   factors = "treatment",
#'   replicate = "bio_rep"
#' )
#'
#' # Multi-factor import
#' dat <- read_pepdiff(
#'   "data.csv",
#'   id = "peptide",
#'   gene = "gene_id",
#'   value = "total_area",
#'   factors = c("treatment", "timepoint"),
#'   replicate = "bio_rep"
#' )
#' }
read_pepdiff <- function(file, id, gene, value, factors, replicate, tech_rep = NULL) {
  # Capture call for reproducibility
  call <- match.call()

  # Validate inputs
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }

  if (!is.character(factors) || length(factors) == 0) {
    stop("'factors' must be a non-empty character vector", call. = FALSE)
  }

  # Read raw data
  raw_data <- readr::read_csv(file, show_col_types = FALSE)

  # Check required columns exist
  required_cols <- c(id, gene, value, factors, replicate)
  if (!is.null(tech_rep)) {
    required_cols <- c(required_cols, tech_rep)
  }
  missing_cols <- setdiff(required_cols, names(raw_data))
  if (length(missing_cols) > 0) {
    stop(
      "Required columns not found in file: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Build standardized data frame
  # Start with required columns
  data <- tibble::tibble(
    peptide = as.character(raw_data[[id]]),
    gene_id = as.character(raw_data[[gene]]),
    bio_rep = as.character(raw_data[[replicate]]),
    value = as.numeric(raw_data[[value]])
  )

  # Add factor columns
  for (f in factors) {
    data[[f]] <- as.character(raw_data[[f]])
  }

  # Add tech_rep if present
  if (!is.null(tech_rep)) {
    data$tech_rep <- as.character(raw_data[[tech_rep]])
  }

  # Remove duplicate rows
  data <- dplyr::distinct(data)

  # Validate factors have multiple levels
  validate_factors(data, factors)

  # Compute design summary
  design <- compute_design_summary(data, factors)

  # Compute missingness
  missingness <- compute_missingness(data)

  # Get unique peptides
  peptides <- unique(data$peptide)

  # Create and validate object
  obj <- new_pepdiff_data(
    data = data,
    factors = factors,
    design = design,
    missingness = missingness,
    peptides = peptides,
    call = call
  )

  validate_pepdiff_data(obj)
}


# =============================================================================
# Technical Replicate Handling
# =============================================================================

#' Combine technical replicates
#'
#' Explicitly combines technical replicates by averaging values within each
#' combination of peptide, factors, and biological replicate.
#'
#' @param data A pepdiff_data object with a tech_rep column
#' @param fun Function to use for combining (default: mean)
#'
#' @return A pepdiff_data object with technical replicates combined
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Import data with technical replicates
#' dat <- read_pepdiff(..., tech_rep = "tech_rep")
#'
#' # Combine by averaging
#' dat <- combine_tech_reps(dat)
#'
#' # Or combine by taking median
#' dat <- combine_tech_reps(dat, fun = median)
#' }
combine_tech_reps <- function(data, fun = mean) {
  UseMethod("combine_tech_reps")
}


#' @export
combine_tech_reps.pepdiff_data <- function(data, fun = mean) {
  if (!"tech_rep" %in% names(data$data)) {
    warning("No 'tech_rep' column found - data returned unchanged", call. = FALSE)
    return(data)
  }

  # Columns to group by (all except tech_rep and value)
  group_cols <- c("peptide", "gene_id", data$factors, "bio_rep")

  # Combine technical replicates
  combined <- data$data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarize(
      value = fun(.data$value, na.rm = TRUE),
      .groups = "drop"
    )

  # Recompute design and missingness
  design <- compute_design_summary(combined, data$factors)
  missingness <- compute_missingness(combined)

  # Create new object
  new_pepdiff_data(
    data = combined,
    factors = data$factors,
    design = design,
    missingness = missingness,
    peptides = data$peptides,
    call = data$call
  )
}


#' @export
combine_tech_reps.data.frame <- function(data, fun = mean) {
  # Legacy support: handle raw data frames
  # Assumes legacy column names: gene_id, peptide, treatment, seconds, bio_rep, quant
  if (!all(c("gene_id", "peptide", "treatment", "seconds", "bio_rep", "quant") %in% names(data))) {
    stop("Data frame must have columns: gene_id, peptide, treatment, seconds, bio_rep, quant",
         call. = FALSE)
  }

  data %>%
    dplyr::group_by(.data$gene_id, .data$peptide, .data$treatment, .data$seconds, .data$bio_rep) %>%
    dplyr::summarize(mean_tr_quant = fun(.data$quant, na.rm = TRUE), .groups = "drop")
}


# =============================================================================
# Print and Summary Methods
# =============================================================================

#' Print method for pepdiff_data
#'
#' @param x A pepdiff_data object
#' @param ... Additional arguments (ignored)
#' @return The object invisibly
#' @export
print.pepdiff_data <- function(x, ...) {
  cat("pepdiff_data object\n")
  cat("-------------------\n")

  # Dimensions
  n_peptides <- length(x$peptides)
  n_obs <- nrow(x$data)
  cat(sprintf("Peptides: %d\n", n_peptides))
  cat(sprintf("Observations: %d\n", n_obs))

  # Factors
  cat(sprintf("Factors: %s\n", paste(x$factors, collapse = ", ")))

  # Design summary
  cat("\nDesign:\n")
  for (i in seq_len(nrow(x$design))) {
    row <- x$design[i, ]
    factor_str <- paste(
      sapply(x$factors, function(f) paste0(f, "=", row[[f]])),
      collapse = ", "
    )
    cat(sprintf("  %s: %d reps\n", factor_str, row$n_reps))
  }

  # Missingness summary
  overall_na_rate <- mean(x$missingness$na_rate) * 100
  if (overall_na_rate > 0) {
    cat(sprintf("\nMissingness: %.1f%% overall\n", overall_na_rate))
  } else {
    cat("\nNo missing values\n")
  }

  invisible(x)
}


#' Summary method for pepdiff_data
#'
#' @param object A pepdiff_data object
#' @param ... Additional arguments (ignored)
#' @return A summary list invisibly
#' @export
summary.pepdiff_data <- function(object, ...) {
  cat("pepdiff_data Summary\n")
  cat("====================\n\n")

  # Basic info
  cat("Data dimensions:\n")
  cat(sprintf("  - Total observations: %d\n", nrow(object$data)))
  cat(sprintf("  - Unique peptides: %d\n", length(object$peptides)))
  cat(sprintf("  - Unique genes: %d\n", dplyr::n_distinct(object$data$gene_id)))

  # Factor breakdown
  cat("\nExperimental factors:\n")
  for (f in object$factors) {
    levels <- unique(object$data[[f]])
    cat(sprintf("  - %s: %s\n", f, paste(levels, collapse = ", ")))
  }

  # Design table
  cat("\nDesign summary:\n")
  print(object$design, n = Inf)

  # Missingness
  cat("\nMissingness:\n")
  n_with_missing <- sum(object$missingness$na_rate > 0)
  cat(sprintf("  - Peptides with missing values: %d (%.1f%%)\n",
              n_with_missing, 100 * n_with_missing / length(object$peptides)))

  if (n_with_missing > 0) {
    cat(sprintf("  - Mean NA rate: %.1f%%\n", 100 * mean(object$missingness$na_rate)))
    cat(sprintf("  - Max NA rate: %.1f%%\n", 100 * max(object$missingness$na_rate)))
  }

  invisible(list(
    n_observations = nrow(object$data),
    n_peptides = length(object$peptides),
    n_genes = dplyr::n_distinct(object$data$gene_id),
    factors = object$factors,
    design = object$design,
    missingness_summary = list(
      n_with_missing = n_with_missing,
      mean_na_rate = mean(object$missingness$na_rate),
      max_na_rate = max(object$missingness$na_rate)
    )
  ))
}


# =============================================================================
# Subset method
# =============================================================================

#' Subset a pepdiff_data object
#'
#' @param x A pepdiff_data object
#' @param peptides Character vector of peptide IDs to keep (optional)
#' @param ... Additional subsetting expressions evaluated in the data
#' @return A new pepdiff_data object
#' @export
subset.pepdiff_data <- function(x, peptides = NULL, ...) {
  data <- x$data

  # Filter by peptide IDs if provided
  if (!is.null(peptides)) {
    data <- data[data$peptide %in% peptides, ]
  }

  # Apply additional filtering via ...
  dots <- rlang::enquos(...)
  if (length(dots) > 0) {
    data <- dplyr::filter(data, !!!dots)
  }

  # Recompute derived components
  design <- compute_design_summary(data, x$factors)
  missingness <- compute_missingness(data)
  new_peptides <- unique(data$peptide)

  new_pepdiff_data(
    data = data,
    factors = x$factors,
    design = design,
    missingness = missingness,
    peptides = new_peptides,
    call = x$call
  )
}
