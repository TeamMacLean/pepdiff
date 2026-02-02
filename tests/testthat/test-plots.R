# Tests for plots.R - visualization functions

# =============================================================================
# Test plot.pepdiff_data
# =============================================================================

test_that("plot.pepdiff_data returns a plot object", {
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

  # Should not error
  p <- plot(imported)
  expect_true(inherits(p, "gg") || inherits(p, "gtable") || inherits(p, "ggplot"))
})


# =============================================================================
# Test plot_pca_simple
# =============================================================================

test_that("plot_pca_simple creates PCA plot", {
  test_data <- make_test_data(n_peptides = 20, n_reps = 3, seed = 42)
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

  p <- plot_pca_simple(imported)
  expect_s3_class(p, "ggplot")
})


test_that("plot_pca_simple handles color_by parameter", {
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

  p <- plot_pca_simple(imported, color_by = "timepoint")
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# Test plot_distributions_simple
# =============================================================================

test_that("plot_distributions_simple creates distribution plot", {
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

  p <- plot_distributions_simple(imported)
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# Test plot_missingness_simple
# =============================================================================

test_that("plot_missingness_simple handles no missing data", {
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

  p <- plot_missingness_simple(imported)
  expect_s3_class(p, "ggplot")
})


test_that("plot_missingness_simple handles missing data", {
  test_data <- make_test_data(na_rate = 0.2, seed = 42)
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

  p <- plot_missingness_simple(imported)
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# Test plot.pepdiff_results
# =============================================================================

test_that("plot.pepdiff_results returns a plot object", {
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

  results <- compare(imported, compare = "treatment", ref = "ctrl", method = "pairwise")

  # Should not error
  p <- plot(results)
  expect_true(inherits(p, "gg") || inherits(p, "gtable") || inherits(p, "ggplot"))
})


# =============================================================================
# Test plot_volcano_new
# =============================================================================

test_that("plot_volcano_new creates volcano plot", {
  test_data <- make_test_data(n_peptides = 20, effect_size = 3.0, seed = 42)
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

  results <- compare(imported, compare = "treatment", ref = "ctrl", method = "pairwise")

  p <- plot_volcano_new(results)
  expect_s3_class(p, "ggplot")
})


test_that("plot_volcano_new handles fc_threshold", {
  test_data <- make_test_data(n_peptides = 20, effect_size = 3.0, seed = 42)
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

  results <- compare(imported, compare = "treatment", ref = "ctrl", method = "pairwise")

  p <- plot_volcano_new(results, fc_threshold = 2)
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# Test plot_pvalue_histogram
# =============================================================================

test_that("plot_pvalue_histogram creates histogram", {
  test_data <- make_test_data(n_peptides = 20, seed = 42)
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

  results <- compare(imported, compare = "treatment", ref = "ctrl", method = "pairwise")

  p <- plot_pvalue_histogram(results)
  expect_s3_class(p, "ggplot")
})


# =============================================================================
# Test plot_fc_distribution_new
# =============================================================================

test_that("plot_fc_distribution_new creates histogram", {
  test_data <- make_test_data(n_peptides = 20, seed = 42)
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

  results <- compare(imported, compare = "treatment", ref = "ctrl", method = "pairwise")

  p <- plot_fc_distribution_new(results)
  expect_s3_class(p, "ggplot")
})
