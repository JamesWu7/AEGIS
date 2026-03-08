test_that("load_10x_lymphnode builds Seurat object from raw files", {
  raw_dir <- normalizePath(file.path(testthat::test_path("..", "..")), mustWork = TRUE)

  req <- c(
    "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5",
    "V1_Human_Lymph_Node_spatial.tar.gz",
    "V1_Human_Lymph_Node_metrics_summary.csv"
  )
  if (!all(file.exists(file.path(raw_dir, req)))) {
    skip("Authoritative raw files are not available in repository root.")
  }

  tmp <- file.path(tempdir(), "aegis_loader_test")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  for (f in req) file.copy(file.path(raw_dir, f), file.path(tmp, f), overwrite = TRUE)

  seu <- load_10x_lymphnode(data_dir = tmp)
  expect_true(inherits(seu, "Seurat"))
  expect_true(ncol(seu) > 0)
  expect_true(dir.exists(file.path(tmp, "spatial")))
})
