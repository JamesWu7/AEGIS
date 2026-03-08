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
