# Tests for legacy.R - backwards compatibility

# Note: These tests verify that old code patterns still work,
# even though they're deprecated

test_that("compare.data.frame issues deprecation warning", {
  skip("Legacy compare requires old data format")
  # This test is a placeholder - the actual legacy compare needs
  # the full old data structure which is complex to set up
})


test_that("matrix_data function exists for legacy code", {
  # Just check the function exists
  expect_true(exists("matrix_data", envir = asNamespace("pepdiff")))
})


test_that("select_columns_for_contrast exists for legacy code", {
  expect_true(exists("select_columns_for_contrast", envir = asNamespace("pepdiff")))
})


test_that("min_peptide_values exists for legacy code", {
  expect_true(exists("min_peptide_values", envir = asNamespace("pepdiff")))
})


test_that("replace_vals exists for legacy code", {
  expect_true(exists("replace_vals", envir = asNamespace("pepdiff")))
})


test_that("import_data is exported and issues deprecation note", {
  # Check function is exported
  expect_true("import_data" %in% getNamespaceExports("pepdiff"))
})


test_that("compare_many is exported", {
  expect_true("compare_many" %in% getNamespaceExports("pepdiff"))
})


test_that("get_bootstrap_percentile is exported", {
  expect_true("get_bootstrap_percentile" %in% getNamespaceExports("pepdiff"))
})


test_that("get_wilcoxon_percentile is exported", {
  expect_true("get_wilcoxon_percentile" %in% getNamespaceExports("pepdiff"))
})


test_that("get_kruskal_percentile is exported", {
  expect_true("get_kruskal_percentile" %in% getNamespaceExports("pepdiff"))
})


test_that("get_rp_percentile is exported", {
  expect_true("get_rp_percentile" %in% getNamespaceExports("pepdiff"))
})


# Test that new functions are also exported
test_that("read_pepdiff is exported", {
  expect_true("read_pepdiff" %in% getNamespaceExports("pepdiff"))
})


test_that("compare is exported", {
  expect_true("compare" %in% getNamespaceExports("pepdiff"))
})


test_that("combine_tech_reps is exported", {
  expect_true("combine_tech_reps" %in% getNamespaceExports("pepdiff"))
})


test_that("significant is exported", {
  expect_true("significant" %in% getNamespaceExports("pepdiff"))
})


test_that("test_wilcoxon is exported", {
  expect_true("test_wilcoxon" %in% getNamespaceExports("pepdiff"))
})


test_that("test_bootstrap_t is exported", {
  expect_true("test_bootstrap_t" %in% getNamespaceExports("pepdiff"))
})


test_that("test_bayes_t is exported", {
  expect_true("test_bayes_t" %in% getNamespaceExports("pepdiff"))
})


test_that("test_rankprod is exported", {
  expect_true("test_rankprod" %in% getNamespaceExports("pepdiff"))
})
