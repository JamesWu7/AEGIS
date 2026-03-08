test_that("compare_methods stores pairwise and spot-level agreement", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 29)

  # Create partial overlap to verify overlap logic.
  deconv$SPOTlight <- deconv$SPOTlight[, setdiff(colnames(deconv$SPOTlight), "Stromal"), drop = FALSE]

  obj <- as_aegis(aegis_example, deconv)

  obj <- compare_methods(obj)
  expect_true("comparison" %in% names(obj$consensus))
  cmp <- obj$consensus$comparison

  expect_true(is.data.frame(cmp$pairwise_celltype_cor))
  expect_true(is.data.frame(cmp$celltype_agreement))
  expect_true(is.data.frame(cmp$spot_agreement))
  expect_true(is.data.frame(cmp$summary))
  expect_true(is.matrix(cmp$heatmap_matrix))

  expect_true(all(c("method_1", "method_2", "celltype", "correlation") %in% colnames(cmp$pairwise_celltype_cor)))
  expect_true(all(c("spot", "agreement") %in% colnames(cmp$spot_agreement)))
  expect_false(any(cmp$pairwise_celltype_cor$method_1 == "SPOTlight" &
    cmp$pairwise_celltype_cor$celltype == "Stromal"))
})

test_that("compute_consensus stores consensus outputs and shared celltypes", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 31)
  deconv$cell2location <- deconv$cell2location[, setdiff(colnames(deconv$cell2location), "NK_cell"), drop = FALSE]

  obj <- as_aegis(aegis_example, deconv)
  obj <- compute_consensus(obj)
  expect_true("result" %in% names(obj$consensus))
  res <- obj$consensus$result

  expect_true(all(c("consensus_matrix", "method_disagreement", "spot_confidence", "celltype_stability", "shared_celltypes") %in% names(res)))
  expect_true(is.matrix(res$consensus_matrix))
  expect_true(is.matrix(res$method_disagreement))
  expect_true(is.data.frame(res$spot_confidence))
  expect_true(is.data.frame(res$celltype_stability))
  expect_true(is.character(res$shared_celltypes))

  expect_false("NK_cell" %in% res$shared_celltypes)
  expect_identical(colnames(res$consensus_matrix), res$shared_celltypes)
  expect_equal(nrow(res$spot_confidence), ncol(aegis_example))
})
