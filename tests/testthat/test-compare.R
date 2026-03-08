test_that("compare_methods and compute_consensus run successfully", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 29)
  obj <- as_aegis(aegis_example, deconv, markers = aegis_default_markers())

  obj <- compare_methods(obj)
  expect_true("comparison" %in% names(obj$consensus))
  expect_true(is.matrix(obj$consensus$comparison$heatmap_matrix))

  obj <- compute_consensus(obj)
  expect_true("result" %in% names(obj$consensus))
  expect_true(all(c("consensus_matrix", "method_disagreement", "spot_confidence", "celltype_stability") %in% names(obj$consensus$result)))
})
