test_that("run_aegis supports single-sample inputs", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 707)

  obj <- run_aegis(
    aegis_example,
    deconv = deconv,
    markers = aegis_default_markers()
  )

  expect_s3_class(obj, "aegis")
  expect_true("basic" %in% names(obj$audit))
  expect_true("comparison" %in% names(obj$consensus))
  expect_true("result" %in% names(obj$consensus))
})

test_that("run_aegis supports multi-sample inputs", {
  data("aegis_example", package = "AEGIS")
  spots <- colnames(aegis_example)
  n_half <- floor(length(spots) / 2)
  seu_list <- list(
    sample1 = suppressWarnings(aegis_example[, spots[seq_len(n_half)]]),
    sample2 = suppressWarnings(aegis_example[, spots[seq.int(n_half + 1L, length(spots))]])
  )

  deconv_nested <- list(
    sample1 = simulate_deconv_results(seu_list$sample1, methods = c("RCTD", "SPOTlight"), seed = 708),
    sample2 = simulate_deconv_results(seu_list$sample2, methods = c("RCTD", "SPOTlight"), seed = 709)
  )

  obj <- run_aegis(seu_list, deconv = deconv_nested, markers = aegis_default_markers())
  expect_s3_class(obj, "aegis_multi")
  expect_true("basic" %in% names(obj$audit))
  expect_true("comparison" %in% names(obj$consensus))
  expect_true("result" %in% names(obj$consensus))
})

test_that("run_aegis accepts an existing aegis object", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 710)
  obj0 <- as_aegis(aegis_example, deconv = deconv)

  obj <- run_aegis(obj0, do_marker = FALSE)
  expect_s3_class(obj, "aegis")
  expect_true("basic" %in% names(obj$audit))
  expect_false("marker" %in% names(obj$audit))
})
