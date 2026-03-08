test_that("audit_basic and audit_marker run and store required fields", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 17)
  obj <- as_aegis(aegis_example, deconv, markers = aegis_default_markers())

  obj <- audit_basic(obj)
  expect_true("basic" %in% names(obj$audit))
  expect_true(all(c("spot_metrics", "summary_table", "sparsity_summary") %in% names(obj$audit$basic)))

  obj <- audit_marker(obj)
  expect_true("marker" %in% names(obj$audit))
  expect_true(all(c("concordance", "summary_table", "marker_usage") %in% names(obj$audit$marker)))

  obj <- audit_spatial(obj)
  expect_true("spatial" %in% names(obj$audit))
  expect_true(all(c("detail", "spot_summary", "method_summary") %in% names(obj$audit$spatial)))
})
