fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

make_multi_fixture <- function() {
  data("aegis_example", package = "AEGIS")
  spots <- colnames(aegis_example)
  n1 <- min(40L, length(spots) %/% 2L)
  n2 <- min(40L, length(spots) - n1)
  s1_cells <- spots[seq_len(n1)]
  s2_cells <- spots[seq.int(from = n1 + 1L, length.out = n2)]

  seu1 <- suppressWarnings(aegis_example[, s1_cells])
  seu2 <- suppressWarnings(aegis_example[, s2_cells])
  seu1 <- AEGIS:::attach_sample_metadata(seu1, sample_id = "sample1", section_id = "section1")
  seu2 <- AEGIS:::attach_sample_metadata(seu2, sample_id = "sample2", section_id = "section2")

  de1 <- simulate_deconv_results(seu1, methods = c("RCTD", "SPOTlight"), seed = 501)
  de2 <- simulate_deconv_results(seu2, methods = c("RCTD", "SPOTlight"), seed = 502)

  list(
    seu_list = list(sample1 = seu1, sample2 = seu2),
    deconv = list(sample1 = de1, sample2 = de2)
  )
}

link_or_copy <- function(src, dst) {
  ok <- suppressWarnings(file.symlink(src, dst))
  if (!isTRUE(ok)) {
    ok <- file.copy(src, dst, overwrite = TRUE)
  }
  ok
}

test_that("load_10x_spatial_set loads multiple sample directories", {
  raw_dir <- normalizePath(file.path(testthat::test_path("..", "..")), mustWork = TRUE)
  req <- c(
    "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5",
    "V1_Human_Lymph_Node_spatial.tar.gz",
    "V1_Human_Lymph_Node_metrics_summary.csv"
  )
  if (!all(file.exists(file.path(raw_dir, req)))) {
    skip("Authoritative raw files are not available in repository root.")
  }

  root <- file.path(tempdir(), "aegis_multi_loader")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  d1 <- file.path(root, "sample1")
  d2 <- file.path(root, "sample2")
  dir.create(d1, recursive = TRUE, showWarnings = FALSE)
  dir.create(d2, recursive = TRUE, showWarnings = FALSE)

  for (f in req) {
    expect_true(link_or_copy(file.path(raw_dir, f), file.path(d1, f)))
    expect_true(link_or_copy(file.path(raw_dir, f), file.path(d2, f)))
  }

  seu_list <- load_10x_spatial_set(paths = c(d1, d2), sample_ids = c("sample1", "sample2"))
  expect_true(is.list(seu_list))
  expect_true(all(vapply(seu_list, inherits, logical(1), what = "Seurat")))
  expect_identical(names(seu_list), c("sample1", "sample2"))
  expect_true(all(seu_list$sample1$sample_id == "sample1"))
  expect_true(all(seu_list$sample2$sample_id == "sample2"))

  expect_error(load_10x_spatial_set(paths = c(d1, file.path(root, "missing"))), "do not exist")
})

test_that("merge_spatial_seurat_list keeps sample traceability and downstream usability", {
  fx <- make_multi_fixture()
  merged <- merge_spatial_seurat_list(fx$seu_list, sample_ids = c("sample1", "sample2"))
  expect_true(inherits(merged, "Seurat"))
  expect_true("sample_id" %in% colnames(merged@meta.data))
  expect_equal(anyDuplicated(colnames(merged)), 0L)

  deconv_merged <- simulate_deconv_results(merged, methods = c("RCTD", "SPOTlight"), seed = 601)
  obj <- as_aegis(merged, deconv = deconv_merged, meta = list(group_var = "sample_id"))
  obj <- audit_basic(obj)
  expect_true(inherits(obj, "aegis"))
  expect_true("by_sample" %in% names(obj$audit$basic))
})

test_that("as_aegis_multi validates nested structures", {
  fx <- make_multi_fixture()
  obj <- as_aegis_multi(fx$seu_list, deconv = fx$deconv, markers = aegis_default_markers())
  expect_true(inherits(obj, "aegis_multi"))
  expect_identical(sort(obj$meta$sample_ids), c("sample1", "sample2"))

  bad_deconv <- fx$deconv
  names(bad_deconv)[1] <- "wrong_sample"
  expect_error(as_aegis_multi(fx$seu_list, bad_deconv), "missing samples")

  bad_nested <- fx$deconv
  bad_nested$sample1$RCTD <- as.matrix(utils::read.csv(fixture_path("multi_sample_deconv_mismatch.csv"), row.names = 1, check.names = FALSE))
  expect_error(as_aegis_multi(fx$seu_list, bad_nested), "missing")
})

test_that("multi-sample audits and comparisons keep sample boundaries", {
  fx <- make_multi_fixture()
  obj <- as_aegis_multi(fx$seu_list, deconv = fx$deconv, markers = aegis_default_markers())

  obj <- audit_basic(obj)
  obj <- audit_marker(obj)
  obj <- audit_spatial(obj)
  obj <- compare_methods(obj)
  obj <- compute_consensus(obj)

  expect_true("by_sample" %in% names(obj$audit$basic))
  expect_true("sample_id" %in% colnames(obj$audit$basic$summary))
  expect_true("sample_id" %in% colnames(obj$audit$marker$concordance))
  expect_true("sample_id" %in% colnames(obj$audit$spatial$spot_metrics))
  expect_true("sample_id" %in% colnames(obj$consensus$comparison$pairwise_summary))
  expect_true("sample_id" %in% colnames(obj$consensus$result$spot_confidence))

  n1 <- names(obj$audit$spatial$by_sample$sample1$neighbors)
  n2 <- names(obj$audit$spatial$by_sample$sample2$neighbors)
  expect_true(all(n1 %in% colnames(fx$seu_list$sample1)))
  expect_true(all(n2 %in% colnames(fx$seu_list$sample2)))
  expect_false(any(n1 %in% colnames(fx$seu_list$sample2)))
})

test_that("split_aegis_by_sample and summarize_by_sample work for multi and single", {
  fx <- make_multi_fixture()
  obj <- as_aegis_multi(fx$seu_list, deconv = fx$deconv, markers = aegis_default_markers())
  obj <- audit_basic(obj)
  obj <- compare_methods(obj)
  obj <- compute_consensus(obj)

  parts <- split_aegis_by_sample(obj)
  expect_identical(names(parts), c("sample1", "sample2"))
  expect_true(all(vapply(parts, inherits, logical(1), what = "aegis")))

  s_tbl <- summarize_by_sample(obj)
  expect_true("sample_id" %in% colnames(s_tbl))
  expect_true("mean_dominance" %in% colnames(s_tbl))
  expect_true(all(c("sample1", "sample2") %in% unique(s_tbl$sample_id)))

  single <- as_aegis(fx$seu_list$sample1, fx$deconv$sample1)
  single <- audit_basic(single)
  single_tbl <- summarize_by_sample(single)
  expect_true("sample_id" %in% colnames(single_tbl))
  expect_equal(length(unique(single_tbl$sample_id)), 1L)
})

test_that("plot functions are safe for multi-sample objects", {
  fx <- make_multi_fixture()
  obj <- as_aegis_multi(fx$seu_list, deconv = fx$deconv, markers = aegis_default_markers())
  obj <- audit_basic(obj)
  obj <- compare_methods(obj)
  obj <- compute_consensus(obj)

  expect_error(plot_audit(obj, type = "dominance"), "provide `sample`")
  expect_error(plot_compare(obj, type = "heatmap"), "provide `sample`")
  expect_error(plot_method_ranking(obj), "provide `sample`")
  expect_error(plot_disagreement_map(obj), "provide `sample`")
  expect_error(plot_consensus_confidence(obj), "provide `sample`")
  expect_s3_class(plot_audit(obj, type = "dominance", sample = "sample1"), "ggplot")
  expect_s3_class(plot_compare(obj, type = "heatmap", sample = "sample1"), "ggplot")
  expect_s3_class(plot_method_ranking(obj, sample = "sample1"), "ggplot")
  p_dis <- plot_disagreement_map(obj, sample = "sample1")
  p_conf <- plot_consensus_confidence(obj, sample = "sample1")
  expect_true(inherits(p_dis, "ggplot") || inherits(p_dis, "patchwork"))
  expect_true(inherits(p_conf, "ggplot") || inherits(p_conf, "patchwork"))
})

test_that("render_report_batch writes per-sample reports and respects overwrite", {
  fx <- make_multi_fixture()
  obj <- as_aegis_multi(fx$seu_list, deconv = fx$deconv, markers = aegis_default_markers())
  obj <- audit_basic(obj)
  obj <- audit_marker(obj)
  obj <- audit_spatial(obj)
  obj <- compare_methods(obj)
  obj <- compute_consensus(obj)

  out_dir <- tempfile("aegis_batch_reports_")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  res <- render_report_batch(obj, output_dir = out_dir, overwrite = FALSE)
  expect_true(file.exists(res$reports[["sample1"]]))
  expect_true(file.exists(res$reports[["sample2"]]))
  expect_true(file.exists(res$summary_file))

  expect_error(render_report_batch(obj, output_dir = out_dir, overwrite = FALSE), "already exists")
  res2 <- render_report_batch(obj, output_dir = out_dir, overwrite = TRUE)
  expect_true(file.exists(res2$summary_file))
})
