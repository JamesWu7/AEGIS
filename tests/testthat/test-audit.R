test_that("audit_basic stores required metrics", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 17)
  obj <- as_aegis(aegis_example, deconv)

  out <- audit_basic(obj)
  expect_s3_class(out, "aegis")
  expect_true("basic" %in% names(out$audit))

  basic <- out$audit$basic
  expect_true(is.list(basic$per_method_spot_metrics))
  expect_true(is.data.frame(basic$summary))

  req_cols <- c("spot", "method", "sum_dev", "dominance", "entropy", "n_detected_types")
  for (m in names(basic$per_method_spot_metrics)) {
    tbl <- basic$per_method_spot_metrics[[m]]
    expect_true(all(req_cols %in% colnames(tbl)))
  }

  expect_true(all(c("method", "zero_fraction", "near_zero_fraction", "mean_dominance", "mean_entropy", "mean_n_detected_types", "mean_sum_dev") %in% colnames(basic$summary)))
})

test_that("audit_basic handles zeros and threshold behavior", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 18)

  # Inject exact zeros while preserving row sums.
  m <- deconv[[1]]
  m[, 1] <- pmax(0, m[, 1] - 0.15)
  rs <- rowSums(m)
  m <- m / rs
  deconv[[1]] <- m

  obj <- as_aegis(aegis_example, deconv)

  low_thr <- audit_basic(obj, threshold = 0.01)$audit$basic
  high_thr <- audit_basic(obj, threshold = 0.30)$audit$basic

  mean_low <- mean(low_thr$spot_metrics$n_detected_types)
  mean_high <- mean(high_thr$spot_metrics$n_detected_types)

  expect_true(is.finite(mean(low_thr$spot_metrics$entropy)))
  expect_gte(mean_low, mean_high)
})

test_that("audit_marker stores concordance and marker usage", {
  data("aegis_example", package = "AEGIS")
  markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))
  deconv <- simulate_deconv_results(aegis_example, seed = 19)
  obj <- as_aegis(aegis_example, deconv, markers = markers)

  out <- audit_marker(obj)
  expect_s3_class(out, "aegis")
  expect_true("marker" %in% names(out$audit))

  marker_res <- out$audit$marker
  expect_true(is.matrix(marker_res$marker_scores))
  expect_true(is.data.frame(marker_res$concordance))
  expect_true(is.data.frame(marker_res$marker_usage))
  expect_true(is.data.frame(marker_res$summary))
  expect_true(all(c("method", "celltype", "pearson_cor", "spearman_cor", "n_markers_used", "n_spots_used") %in% colnames(marker_res$concordance)))
  expect_true(all(c("celltype", "n_markers_requested", "n_markers_used") %in% colnames(marker_res$marker_usage)))
})

test_that("audit_marker handles missing markers/genes gracefully", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 20)
  obj <- as_aegis(aegis_example, deconv, markers = NULL)

  expect_warning(out <- audit_marker(obj), "No marker list")
  expect_true("marker" %in% names(out$audit))
  expect_equal(nrow(out$audit$marker$concordance), 0)

  partial_markers <- list(
    B_cell = c("CD79A", "NOT_A_GENE_1"),
    Unknown = c("NOT_A_GENE_2")
  )
  obj2 <- as_aegis(aegis_example, deconv, markers = partial_markers)
  expect_warning(out2 <- audit_marker(obj2), "no genes present")
  usage <- out2$audit$marker$marker_usage
  expect_true("Unknown" %in% usage$celltype)
  expect_equal(usage$n_markers_used[usage$celltype == "Unknown"], 0)
})

test_that("audit_spatial stores spatial metrics and validates k", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 21)
  obj <- as_aegis(aegis_example, deconv)

  out <- audit_spatial(obj, k = 6)
  expect_s3_class(out, "aegis")
  expect_true("spatial" %in% names(out$audit))

  spatial <- out$audit$spatial
  expect_true(is.data.frame(spatial$spot_metrics))
  expect_true(is.data.frame(spatial$summary))
  expect_true(is.list(spatial$neighbors))
  expect_true(all(c("spot", "method", "local_inconsistency", "smoothness") %in% colnames(spatial$spot_metrics)))
  expect_true(all(c("method", "mean_local_inconsistency", "mean_smoothness") %in% colnames(spatial$summary)))

  expect_error(audit_spatial(obj, k = 0), "positive integer")
  expect_error(audit_spatial(obj, k = ncol(aegis_example)), "must be <=")
})
