# Tests for GLM/ART stratification via `within` parameter
# These tests verify that the `within` parameter correctly produces
# separate results per stratum, not a diluted main effect.

# =============================================================================
# Test Fixtures
# =============================================================================

#' Generate factorial test data with interaction pattern
#'
#' Creates data where some peptides have an effect only at 24h (not at 0h),
#' simulating a time-dependent treatment response (interaction pattern).
#'
#' @param n_peptides Total number of peptides
#' @param n_interaction Number of peptides with interaction pattern (effect at 24h only)
#' @param n_reps Number of biological replicates per condition
#' @param effect_size Fold change for effect at 24h
#' @param seed Random seed
#' @return A tibble suitable for read_pepdiff()
make_interaction_data <- function(n_peptides = 30, n_interaction = 10, n_reps = 4,
                                   effect_size = 4.0, seed = 42) {
  set.seed(seed)

  # Create peptide and gene IDs
  peptide_ids <- paste0("PEP_", sprintf("%03d", seq_len(n_peptides)))
  gene_ids <- paste0("GENE_", sprintf("%03d", ceiling(seq_len(n_peptides) / 2)))

  # Interaction peptides are last n_interaction
  interaction_peptides <- tail(peptide_ids, n_interaction)

  # Design: treatment (ctrl, trt) x timepoint (0h, 24h)
  factor_grid <- expand.grid(
    treatment = c("ctrl", "trt"),
    timepoint = c("0h", "24h"),
    stringsAsFactors = FALSE
  )

  base_mean <- 1e6
  cv <- 0.25

  data_list <- list()

  for (i in seq_len(nrow(factor_grid))) {
    condition <- factor_grid[i, ]

    for (rep_idx in seq_len(n_reps)) {
      for (pep_idx in seq_len(n_peptides)) {
        peptide <- peptide_ids[pep_idx]
        gene <- gene_ids[pep_idx]

        # Determine fold change
        fc_multiplier <- 1.0
        if (peptide %in% interaction_peptides) {
          # Interaction pattern: effect ONLY at 24h, not at 0h
          if (condition$treatment == "trt" && condition$timepoint == "24h") {
            fc_multiplier <- effect_size
          }
        }

        # Generate Gamma-distributed value
        shape <- 1 / (cv^2)
        adjusted_mean <- base_mean * fc_multiplier
        rate <- shape / adjusted_mean
        value <- stats::rgamma(1, shape = shape, rate = rate)

        row <- data.frame(
          peptide = peptide,
          gene_id = gene,
          bio_rep = rep_idx,
          value = value,
          treatment = condition$treatment,
          timepoint = condition$timepoint,
          stringsAsFactors = FALSE
        )

        data_list[[length(data_list) + 1]] <- row
      }
    }
  }

  tibble::as_tibble(do.call(rbind, data_list))
}


# =============================================================================
# Test 1: GLM with within produces results for each stratum
# =============================================================================

test_that("GLM with within produces results for each stratum", {
  # Setup: 2x2 factorial data
  raw_data <- make_factorial_test_data(seed = 111)

  # Write and read to create pepdiff_data
  temp_file <- write_test_csv(raw_data)
  dat <- read_pepdiff(
    file = temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint"),
    replicate = "bio_rep"
  )
  unlink(temp_file)

  # Run stratified analysis
  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "glm")

  # Should have 2 rows per peptide (one per timepoint)
  expect_equal(nrow(results$results), 2 * length(dat$peptides))

  # Should have timepoint column
  expect_true("timepoint" %in% names(results$results))

  # Should have both timepoint levels
  expect_setequal(unique(results$results$timepoint), c("0h", "24h"))
})


# =============================================================================
# Test 2: Stratified results detect interaction correctly
# =============================================================================

test_that("GLM stratification detects interaction pattern", {
  # Setup: data with interaction pattern (Group C: effect at 24h only)
  raw_data <- make_interaction_data(
    n_peptides = 30,
    n_interaction = 10,  # PEP_021 through PEP_030 have interaction pattern
    n_reps = 4,
    effect_size = 4.0,
    seed = 222
  )

  # Write and read to create pepdiff_data
  temp_file <- write_test_csv(raw_data)
  dat <- read_pepdiff(
    file = temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint"),
    replicate = "bio_rep"
  )
  unlink(temp_file)

  # Run stratified analysis
  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "glm")

  # Interaction peptides (effect at 24h only)
  interaction_peps <- paste0("PEP_", sprintf("%03d", 21:30))

  interaction_results <- results$results |>
    dplyr::filter(peptide %in% interaction_peps) |>
    dplyr::group_by(timepoint) |>
    dplyr::summarise(median_fc = median(fold_change), .groups = "drop")

  # At 0h: FC should be ~1 (no effect)
  fc_0h <- interaction_results$median_fc[interaction_results$timepoint == "0h"]
  expect_true(fc_0h > 0.7 && fc_0h < 1.5,
              label = sprintf("FC at 0h should be ~1, got %.2f", fc_0h))

  # At 24h: FC should be ~4 (the true effect)
  fc_24h <- interaction_results$median_fc[interaction_results$timepoint == "24h"]
  expect_true(fc_24h > 3 && fc_24h < 5,
              label = sprintf("FC at 24h should be ~4, got %.2f", fc_24h))
})


# =============================================================================
# Test 3: FDR correction is applied within each stratum
# =============================================================================

test_that("FDR correction is applied within each stratum", {
  # Setup
  raw_data <- make_factorial_test_data(seed = 333)

  temp_file <- write_test_csv(raw_data)
  dat <- read_pepdiff(
    file = temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint"),
    replicate = "bio_rep"
  )
  unlink(temp_file)

  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "glm")

  # Check FDR is computed separately per stratum
  for (tp in c("0h", "24h")) {
    stratum_results <- results$results |> dplyr::filter(timepoint == tp)
    expected_fdr <- p.adjust(stratum_results$p_value, method = "BH")
    expect_equal(stratum_results$fdr, expected_fdr,
                 label = paste("FDR within", tp))
  }
})


# =============================================================================
# Test 4: ART method with stratification
# =============================================================================

test_that("ART with within produces results for each stratum", {
  skip_if_not_installed("ARTool")

  # Setup
  raw_data <- make_factorial_test_data(seed = 444)

  temp_file <- write_test_csv(raw_data)
  dat <- read_pepdiff(
    file = temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint"),
    replicate = "bio_rep"
  )
  unlink(temp_file)

  # Run stratified ART analysis
  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "art")

  # Should have 2 rows per peptide (one per timepoint)
  expect_equal(nrow(results$results), 2 * length(dat$peptides))

  # Should have timepoint column
  expect_true("timepoint" %in% names(results$results))

  # Should have both timepoint levels
  expect_setequal(unique(results$results$timepoint), c("0h", "24h"))
})


# =============================================================================
# Test 5: Pairwise method with stratification still works
# =============================================================================

test_that("Pairwise method with within still works", {
  raw_data <- make_factorial_test_data(seed = 555)

  temp_file <- write_test_csv(raw_data)
  dat <- read_pepdiff(
    file = temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint"),
    replicate = "bio_rep"
  )
  unlink(temp_file)

  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "pairwise", test = "wilcoxon")

  # Should have timepoint column
  expect_true("timepoint" %in% names(results$results))

  # Should have 2 rows per peptide
  expect_equal(nrow(results$results), 2 * length(dat$peptides))

  # Should have both timepoint levels
  expect_setequal(unique(results$results$timepoint), c("0h", "24h"))
})


# =============================================================================
# Test 6: Pairwise bayes_t with stratification works (bf_threshold passed)
# =============================================================================

test_that("Pairwise bayes_t with within works", {
  raw_data <- make_factorial_test_data(seed = 666)

  temp_file <- write_test_csv(raw_data)
  dat <- read_pepdiff(
    file = temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint"),
    replicate = "bio_rep"
  )
  unlink(temp_file)

  # Run with bayes_t - bf_threshold should be passed through
  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = "timepoint", method = "pairwise",
                     test = "bayes_t", bf_threshold = 10)

  # Should have timepoint column
  expect_true("timepoint" %in% names(results$results))

  # Should have bf column (bayes_t specific)
  expect_true("bf" %in% names(results$results))

  # Should have both timepoint levels
  expect_setequal(unique(results$results$timepoint), c("0h", "24h"))
})


# =============================================================================
# Test 7: Multiple within factors
# =============================================================================

test_that("Multiple within factors work", {
  # Generate 3-factor data: treatment x timepoint x tissue
  set.seed(777)
  factor_grid <- expand.grid(
    treatment = c("ctrl", "trt"),
    timepoint = c("0h", "24h"),
    tissue = c("liver", "kidney"),
    stringsAsFactors = FALSE
  )

  n_peptides <- 5
  n_reps <- 2
  base_mean <- 1e6
  cv <- 0.3

  peptide_ids <- paste0("PEP_", sprintf("%03d", seq_len(n_peptides)))
  gene_ids <- paste0("GENE_", sprintf("%03d", ceiling(seq_len(n_peptides) / 2)))

  data_list <- list()
  for (i in seq_len(nrow(factor_grid))) {
    condition <- factor_grid[i, ]
    for (rep_idx in seq_len(n_reps)) {
      for (pep_idx in seq_len(n_peptides)) {
        shape <- 1 / (cv^2)
        rate <- shape / base_mean
        value <- stats::rgamma(1, shape = shape, rate = rate)

        row <- data.frame(
          peptide = peptide_ids[pep_idx],
          gene_id = gene_ids[pep_idx],
          bio_rep = rep_idx,
          value = value,
          treatment = condition$treatment,
          timepoint = condition$timepoint,
          tissue = condition$tissue,
          stringsAsFactors = FALSE
        )
        data_list[[length(data_list) + 1]] <- row
      }
    }
  }
  raw_data <- tibble::as_tibble(do.call(rbind, data_list))

  temp_file <- write_test_csv(raw_data)
  dat <- read_pepdiff(
    file = temp_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint", "tissue"),
    replicate = "bio_rep"
  )
  unlink(temp_file)

  # Stratify by both timepoint AND tissue
  results <- compare(dat, compare = "treatment", ref = "ctrl",
                     within = c("timepoint", "tissue"), method = "glm")

  # Should have 4 rows per peptide (2 timepoints x 2 tissues)
  expect_equal(nrow(results$results), 4 * length(dat$peptides))

  # Should have both within columns
  expect_true("timepoint" %in% names(results$results))
  expect_true("tissue" %in% names(results$results))

  # Should have all 4 combinations
  combos <- unique(results$results[, c("timepoint", "tissue")])
  expect_equal(nrow(combos), 4)
})
