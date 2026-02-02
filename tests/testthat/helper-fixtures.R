# Test fixture generators for pepdiff
# These functions create synthetic proteomics data for testing

#' Generate synthetic proteomics test data
#'
#' Creates realistic PRM proteomics data with configurable parameters.
#' Data follows Gamma distribution typical of mass spec abundances.
#'
#' @param n_peptides Number of peptides to generate
#' @param n_reps Number of biological replicates per condition
#' @param factors Named list of factor levels, e.g., list(treatment = c("ctrl", "trt"))
#' @param effect_size Fold change for affected peptides
#' @param affected_fraction Proportion of peptides with differential abundance
#' @param na_rate Proportion of missing values to introduce
#' @param seed Random seed for reproducibility
#' @return A tibble suitable for read_pepdiff() or testing
make_test_data <- function(
    n_peptides = 10,
    n_reps = 3,
    factors = list(treatment = c("ctrl", "trt")),
    effect_size = 2.0,
    affected_fraction = 0.3,
    na_rate = 0.05,
    seed = 42
) {
  set.seed(seed)


  # Generate all factor combinations
  factor_grid <- expand.grid(factors, stringsAsFactors = FALSE)

  # Create peptide and gene IDs
  peptide_ids <- paste0("PEP_", sprintf("%03d", seq_len(n_peptides)))
  gene_ids <- paste0("GENE_", sprintf("%03d", ceiling(seq_len(n_peptides) / 2)))

  # Determine which peptides are affected (have differential abundance)
  n_affected <- round(n_peptides * affected_fraction)
  affected_peptides <- sample(peptide_ids, n_affected)

  # Base abundance parameters (realistic proteomics scale)
  base_mean <- 1e6
  cv <- 0.3  # Coefficient of variation

  # Build the full dataset
  data_list <- list()

  for (i in seq_len(nrow(factor_grid))) {
    condition <- factor_grid[i, , drop = FALSE]

    for (rep_idx in seq_len(n_reps)) {
      for (pep_idx in seq_len(n_peptides)) {
        peptide <- peptide_ids[pep_idx]
        gene <- gene_ids[pep_idx]

        # Determine if this peptide is affected and if this is treatment
        is_affected <- peptide %in% affected_peptides

        # Calculate fold change multiplier
        # Apply effect to non-reference levels of the first factor
        fc_multiplier <- 1.0
        if (is_affected) {
          first_factor <- names(factors)[1]
          ref_level <- factors[[first_factor]][1]
          current_level <- condition[[first_factor]]
          if (current_level != ref_level) {
            fc_multiplier <- effect_size
          }
        }

        # Generate Gamma-distributed value
        # shape = 1/CV^2, rate = shape/mean
        shape <- 1 / (cv^2)
        adjusted_mean <- base_mean * fc_multiplier
        rate <- shape / adjusted_mean
        value <- stats::rgamma(1, shape = shape, rate = rate)

        # Create row
        row <- data.frame(
          peptide = peptide,
          gene_id = gene,
          bio_rep = rep_idx,
          value = value,
          stringsAsFactors = FALSE
        )

        # Add factor columns
        for (f in names(factors)) {
          row[[f]] <- condition[[f]]
        }

        data_list[[length(data_list) + 1]] <- row
      }
    }
  }

  result <- do.call(rbind, data_list)

  # Introduce missing values (MCAR pattern for simplicity)
  if (na_rate > 0) {
    n_missing <- round(nrow(result) * na_rate)
    if (n_missing > 0) {
      missing_idx <- sample(seq_len(nrow(result)), n_missing)
      result$value[missing_idx] <- NA
    }
  }

  tibble::as_tibble(result)
}


#' Generate minimal test data for quick tests
#'
#' Creates a small dataset with 5 peptides and 2 reps per condition
#'
#' @param seed Random seed
#' @return A tibble
make_minimal_test_data <- function(seed = 123) {
  make_test_data(
    n_peptides = 5,
    n_reps = 2,
    factors = list(treatment = c("ctrl", "trt")),
    effect_size = 2.0,
    affected_fraction = 0.4,
    na_rate = 0,
    seed = seed
  )
}


#' Generate multi-factor test data
#'
#' Creates data with two factors (treatment x timepoint) for factorial design tests
#'
#' @param seed Random seed
#' @return A tibble
make_factorial_test_data <- function(seed = 456) {
  make_test_data(
    n_peptides = 10,
    n_reps = 3,
    factors = list(
      treatment = c("ctrl", "trt"),
      timepoint = c("0h", "24h")
    ),
    effect_size = 2.0,
    affected_fraction = 0.3,
    na_rate = 0.05,
    seed = seed
  )
}


#' Generate test data with technical replicates
#'
#' @param n_tech_reps Number of technical replicates per biological replicate
#' @param seed Random seed
#' @return A tibble with tech_rep column
make_tech_rep_test_data <- function(n_tech_reps = 2, seed = 789) {
  base_data <- make_test_data(
    n_peptides = 5,
    n_reps = 2,
    factors = list(treatment = c("ctrl", "trt")),
    effect_size = 2.0,
    affected_fraction = 0.4,
    na_rate = 0,
    seed = seed
  )

  # Expand with technical replicates
  set.seed(seed + 1)
  expanded <- list()

  for (i in seq_len(nrow(base_data))) {
    for (tech in seq_len(n_tech_reps)) {
      row <- base_data[i, ]
      # Add small technical variation (CV ~10%)
      row$value <- row$value * stats::rlnorm(1, 0, 0.1)
      row$tech_rep <- tech
      expanded[[length(expanded) + 1]] <- row
    }
  }

  tibble::as_tibble(do.call(rbind, expanded))
}


#' Write test data to a temporary CSV file
#'
#' @param data A tibble of test data
#' @return Path to the temporary file
write_test_csv <- function(data) {
  temp_file <- tempfile(fileext = ".csv")
  readr::write_csv(data, temp_file)
  temp_file
}


#' Generate test data with zeros (for testing Gamma GLM error handling)
#'
#' @param seed Random seed
#' @return A tibble with some zero values
make_zero_test_data <- function(seed = 999) {
  data <- make_minimal_test_data(seed = seed)
  # Replace some values with zeros
  data$value[c(1, 5, 10)] <- 0
  data
}
