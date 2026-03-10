test_that("single-sample simulated workflow runs end-to-end", {
  data("aegis_example", package = "AEGIS")
  markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))
  deconv <- simulate_deconv_results(aegis_example, seed = 9001)

  obj <- run_aegis(aegis_example, deconv = deconv, markers = markers)
  expect_s3_class(obj, "aegis")
  expect_true(!is.null(obj$audit$basic))
  expect_true(!is.null(obj$audit$marker))
  expect_true(!is.null(obj$audit$spatial))
  expect_true(!is.null(obj$consensus$comparison))
  expect_true(!is.null(obj$consensus$result))

  p1 <- plot_audit(obj, "dominance")
  p2 <- plot_compare(obj, "heatmap")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")

  out <- tempfile(fileext = ".html")
  report <- render_report(obj, output_file = out)
  expect_true(file.exists(report))
  expect_gt(file.info(report)$size, 0)
})

test_that("single-sample imported workflow runs end-to-end", {
  data("aegis_example", package = "AEGIS")
  seu <- suppressWarnings(aegis_example[, colnames(aegis_example)[1:8]])
  spots <- colnames(seu)
  markers <- aegis_default_markers()

  tmp_rctd <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      barcode = spots,
      B_cell = c(0.4, 0.2, 0.1, 0.3, 0.6, 0.1, 0.2, 0.5),
      T_cell = c(0.4, 0.5, 0.3, 0.3, 0.3, 0.7, 0.6, 0.3),
      Myeloid = c(0.2, 0.3, 0.6, 0.4, 0.1, 0.2, 0.2, 0.2),
      check.names = FALSE
    ),
    tmp_rctd,
    row.names = FALSE
  )

  tmp_sp <- tempfile(fileext = ".tsv")
  utils::write.table(
    data.frame(
      spot_id = spots,
      B_cell = c(0.35, 0.2, 0.1, 0.2, 0.5, 0.1, 0.2, 0.45),
      T_cell = c(0.45, 0.5, 0.3, 0.4, 0.3, 0.7, 0.6, 0.35),
      Myeloid = c(0.2, 0.3, 0.6, 0.4, 0.2, 0.2, 0.2, 0.2),
      sample = "S1",
      check.names = FALSE
    ),
    tmp_sp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  tmp_c2l <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      spot = spots,
      B_cell = c(15, 5, 3, 7, 20, 2, 3, 11),
      T_cell = c(14, 12, 8, 8, 11, 17, 12, 10),
      Myeloid = c(4, 6, 9, 9, 3, 5, 5, 6),
      x = seq_along(spots),
      y = rev(seq_along(spots)),
      check.names = FALSE
    ),
    tmp_c2l,
    row.names = FALSE
  )

  obj <- as_aegis(
    seu,
    deconv = list(
      RCTD = read_rctd(tmp_rctd),
      SPOTlight = read_spotlight(tmp_sp),
      cell2location = read_cell2location(tmp_c2l)
    ),
    markers = markers
  )
  obj <- audit_basic(obj)
  obj <- audit_marker(obj)
  obj <- audit_spatial(obj)
  obj <- compare_methods(obj)
  obj <- compute_consensus(obj)

  expect_s3_class(obj, "aegis")
  expect_true(!is.null(obj$audit$basic$summary))
  expect_true(!is.null(obj$consensus$result$consensus_matrix))

  expect_s3_class(plot_compare(obj, "consensus_map"), "ggplot")
})

test_that("multi-sample workflow runs end-to-end with sample boundaries preserved", {
  data("aegis_example", package = "AEGIS")
  markers <- aegis_default_markers()
  spots <- colnames(aegis_example)
  n_half <- floor(length(spots) / 2)

  seu_list <- list(
    sample1 = suppressWarnings(aegis_example[, spots[seq_len(n_half)]]),
    sample2 = suppressWarnings(aegis_example[, spots[seq.int(n_half + 1L, length(spots))]])
  )

  deconv_nested <- list(
    sample1 = simulate_deconv_results(seu_list$sample1, methods = c("RCTD", "SPOTlight"), seed = 9002),
    sample2 = simulate_deconv_results(seu_list$sample2, methods = c("RCTD", "SPOTlight"), seed = 9003)
  )

  obj_multi <- run_aegis(seu_list, deconv = deconv_nested, markers = markers)
  expect_s3_class(obj_multi, "aegis_multi")

  sum_tbl <- summarize_by_sample(obj_multi)
  expect_true(is.data.frame(sum_tbl))
  expect_true("sample_id" %in% colnames(sum_tbl))
  expect_true(all(c("sample1", "sample2") %in% unique(sum_tbl$sample_id)))

  # Spatial neighborhoods must remain within each sample.
  s1_neighbors <- names(obj_multi$audit$spatial$by_sample$sample1$neighbors)
  s2_neighbors <- names(obj_multi$audit$spatial$by_sample$sample2$neighbors)
  expect_false(any(s1_neighbors %in% colnames(seu_list$sample2)))
  expect_false(any(s2_neighbors %in% colnames(seu_list$sample1)))

  expect_error(plot_audit(obj_multi, "dominance"), "provide `sample`")
  expect_s3_class(plot_audit(obj_multi, "dominance", sample = "sample1"), "ggplot")

  out_dir <- tempfile("aegis_batch_report_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  batch <- render_report_batch(obj_multi, output_dir = out_dir, overwrite = TRUE)
  expect_true(file.exists(batch$reports[["sample1"]]))
  expect_true(file.exists(batch$reports[["sample2"]]))
  expect_true(file.exists(batch$summary_file))
})

test_that("loader + run_aegis full path works when authoritative raw files are present", {
  raw_dir <- normalizePath(file.path(testthat::test_path("..", "..")), mustWork = TRUE)
  req <- c(
    "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5",
    "V1_Human_Lymph_Node_spatial.tar.gz",
    "V1_Human_Lymph_Node_metrics_summary.csv"
  )
  if (!all(file.exists(file.path(raw_dir, req)))) {
    skip("Authoritative raw files are not available in repository root.")
  }

  markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))
  seu <- load_10x_lymphnode(data_dir = raw_dir)
  deconv <- simulate_deconv_results(seu, methods = c("RCTD", "SPOTlight"), seed = 9004)
  obj <- run_aegis(seu, deconv = deconv, markers = markers, do_consensus = TRUE)
  expect_s3_class(obj, "aegis")
  expect_true(!is.null(obj$audit$basic))
})
