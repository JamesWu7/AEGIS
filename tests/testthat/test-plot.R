make_plot_ready_obj <- function() {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 321)
  obj <- as_aegis(aegis_example, deconv, markers = aegis_default_markers())
  obj <- audit_basic(obj)
  obj <- audit_marker(obj)
  obj <- audit_spatial(obj)
  obj <- compare_methods(obj)
  obj <- compute_consensus(obj)
  obj
}

test_that("plot_audit returns ggplot objects for required types", {
  obj <- make_plot_ready_obj()

  expect_s3_class(plot_audit(obj, type = "sumdev"), "ggplot")
  expect_s3_class(plot_audit(obj, type = "dominance", method = "RCTD", palette = "viridis"), "ggplot")
  expect_s3_class(plot_audit(obj, type = "entropy", palette = "brewer"), "ggplot")
  expect_s3_class(plot_audit(obj, type = "marker", palette = "scico"), "ggplot")
  expect_s3_class(plot_audit(obj, type = "smoothness", method = "SPOTlight"), "ggplot")
})

test_that("plot_audit validates missing prerequisites and arguments", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 11)
  obj <- as_aegis(aegis_example, deconv)

  expect_error(plot_audit(obj, type = "sumdev"), "Run audit_basic")
  expect_error(plot_audit(obj, type = "marker"), "Run audit_marker")
  expect_error(plot_audit(obj, type = "smoothness"), "Run audit_spatial")
  expect_error(plot_audit(obj, type = "bad_type"), "should be one of")
})

test_that("plot_compare returns ggplot objects for required types", {
  obj <- make_plot_ready_obj()

  expect_s3_class(plot_compare(obj, type = "heatmap"), "ggplot")
  expect_s3_class(plot_compare(obj, type = "spot_agreement", palette = "viridis"), "ggplot")
  expect_s3_class(plot_compare(obj, type = "consensus_map", palette = "brewer"), "ggplot")
})

test_that("plot_compare validates missing prerequisites and arguments", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 15)
  obj <- as_aegis(aegis_example, deconv)

  expect_error(plot_compare(obj, type = "heatmap"), "Run compare_methods")
  expect_error(plot_compare(obj, type = "spot_agreement"), "Run compare_methods")
  expect_error(plot_compare(obj, type = "consensus_map"), "Run compute_consensus")
  expect_error(plot_compare(obj, type = "invalid"), "should be one of")
})
