# Tests for data.R - pepdiff_data class and import functions

# =============================================================================
# Test pepdiff_data constructor and validator
# =============================================================================

test_that("new_pepdiff_data creates object with correct structure", {
  data <- tibble::tibble(
    peptide = c("P1", "P1", "P2", "P2"),
    gene_id = c("G1", "G1", "G1", "G1"),
    treatment = c("ctrl", "trt", "ctrl", "trt"),
    bio_rep = c("1", "1", "1", "1"),
    value = c(100, 200, 150, 300)
  )

  obj <- new_pepdiff_data(
    data = data,
    factors = "treatment",
    design = tibble::tibble(treatment = c("ctrl", "trt"), n_reps = c(1, 1), n_peptides = c(2, 2), n_observations = c(2, 2)),
    missingness = tibble::tibble(peptide = c("P1", "P2"), na_rate = c(0, 0), mnar_score = c(0, 0), mean_abundance = c(150, 225)),
    peptides = c("P1", "P2"),
    call = quote(test())
  )

  expect_s3_class(obj, "pepdiff_data")
  expect_equal(obj$factors, "treatment")
  expect_equal(length(obj$peptides), 2)
})


test_that("validate_pepdiff_data catches missing components", {
  # Missing required component
  bad_obj <- structure(
    list(data = tibble::tibble(), factors = "a"),
    class = "pepdiff_data"
  )

  expect_error(
    validate_pepdiff_data(bad_obj),
    "missing required components"
  )
})


test_that("validate_pepdiff_data catches missing columns in data", {
  bad_obj <- structure(
    list(
      data = tibble::tibble(peptide = "P1"),  # missing other columns
      factors = "treatment",
      design = tibble::tibble(),
      missingness = tibble::tibble(),
      peptides = "P1",
      call = quote(test())
    ),
    class = "pepdiff_data"
  )

  expect_error(
    validate_pepdiff_data(bad_obj),
    "missing required columns"
  )
})


# =============================================================================
# Test read_pepdiff
# =============================================================================

test_that("read_pepdiff imports CSV correctly", {
  # Create test data and write to temp file
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  result <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  expect_s3_class(result, "pepdiff_data")
  expect_equal(result$factors, "treatment")
  expect_equal(length(result$peptides), 5)
})


test_that("read_pepdiff handles multiple factors", {
  test_data <- make_factorial_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  result <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint"),
    replicate = "bio_rep"
  )

  expect_equal(result$factors, c("treatment", "timepoint"))
  expect_equal(nrow(result$design), 4)  # 2x2 factorial
})


test_that("read_pepdiff errors on missing file", {
  expect_error(
    read_pepdiff(
      file = "nonexistent.csv",
      id = "peptide",
      gene = "gene_id",
      value = "value",
      factors = "treatment",
      replicate = "bio_rep"
    ),
    "File not found"
  )
})


test_that("read_pepdiff errors on missing columns", {
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  expect_error(
    read_pepdiff(
      file = test_file,
      id = "nonexistent_col",
      gene = "gene_id",
      value = "value",
      factors = "treatment",
      replicate = "bio_rep"
    ),
    "not found in file"
  )
})


test_that("read_pepdiff preserves tech_rep column when specified", {
  test_data <- make_tech_rep_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  result <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep",
    tech_rep = "tech_rep"
  )

  expect_true("tech_rep" %in% names(result$data))
})


test_that("read_pepdiff removes duplicate rows", {
  # Create data with duplicates
  test_data <- make_minimal_test_data()
  test_data <- rbind(test_data, test_data[1:2, ])
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  result <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  # Should have removed duplicates
  expect_equal(nrow(result$data), nrow(make_minimal_test_data()))
})


# =============================================================================
# Test combine_tech_reps
# =============================================================================

test_that("combine_tech_reps.pepdiff_data combines technical replicates", {
  test_data <- make_tech_rep_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep",
    tech_rep = "tech_rep"
  )

  # Count rows before combining
  n_before <- nrow(imported$data)

  combined <- combine_tech_reps(imported)

  # Should have fewer rows after combining
  expect_lt(nrow(combined$data), n_before)
  expect_false("tech_rep" %in% names(combined$data))
})


test_that("combine_tech_reps.pepdiff_data warns when no tech_rep column", {
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  expect_warning(
    combine_tech_reps(imported),
    "No 'tech_rep' column found"
  )
})


test_that("combine_tech_reps uses custom function", {
  test_data <- make_tech_rep_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep",
    tech_rep = "tech_rep"
  )

  combined_mean <- combine_tech_reps(imported, fun = mean)
  combined_median <- combine_tech_reps(imported, fun = median)

  # Values should potentially differ with different functions
  expect_s3_class(combined_mean, "pepdiff_data")
  expect_s3_class(combined_median, "pepdiff_data")
})


# =============================================================================
# Test print method
# =============================================================================

test_that("print.pepdiff_data produces output", {
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  output <- capture.output(print(imported))

  expect_true(any(grepl("pepdiff_data", output)))
  expect_true(any(grepl("Peptides:", output)))
  expect_true(any(grepl("Factors:", output)))
})


# =============================================================================
# Test summary method
# =============================================================================

test_that("summary.pepdiff_data produces detailed output", {
  test_data <- make_factorial_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = c("treatment", "timepoint"),
    replicate = "bio_rep"
  )

  output <- capture.output(result <- summary(imported))

  expect_true(any(grepl("Summary", output)))
  expect_true(any(grepl("Experimental factors", output)))
  expect_type(result, "list")
  expect_true("n_peptides" %in% names(result))
})


# =============================================================================
# Test subset method
# =============================================================================

test_that("subset.pepdiff_data filters by peptide IDs", {
  test_data <- make_minimal_test_data()
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  imported <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  # Get first two peptide IDs
  keep_peptides <- imported$peptides[1:2]

  subsetted <- subset(imported, peptides = keep_peptides)

  expect_equal(length(subsetted$peptides), 2)
  expect_true(all(subsetted$data$peptide %in% keep_peptides))
})


# =============================================================================
# Test with missing values
# =============================================================================

test_that("read_pepdiff handles missing values correctly", {
  test_data <- make_test_data(na_rate = 0.1)
  test_file <- write_test_csv(test_data)
  on.exit(unlink(test_file))

  result <- read_pepdiff(
    file = test_file,
    id = "peptide",
    gene = "gene_id",
    value = "value",
    factors = "treatment",
    replicate = "bio_rep"
  )

  expect_s3_class(result, "pepdiff_data")
  expect_true(any(result$missingness$na_rate > 0))
})
