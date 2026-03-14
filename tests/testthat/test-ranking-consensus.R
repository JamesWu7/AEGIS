make_ranking_ready_obj <- function(seed = 901) {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = seed)
  obj <- as_aegis(aegis_example, deconv = deconv, markers = aegis_default_markers())
  obj <- audit_basic(obj)
  obj <- audit_marker(obj)
  obj <- audit_spatial(obj)
  obj <- compare_methods(obj)
  obj
}

test_that("score_methods creates method evidence table with required columns", {
  obj <- make_ranking_ready_obj(seed = 902)
  obj <- score_methods(obj)

  ev <- obj$consensus$method_evidence
  expect_true(is.data.frame(ev))
  expect_true(all(c(
    "method", "marker_score", "spatial_score", "agreement_score",
    "stability_score", "prior_score", "overall_score_raw"
  ) %in% colnames(ev)))
  expect_true(is.numeric(ev$overall_score_raw))
  expect_equal(nrow(ev), length(obj$deconv))
})

test_that("score_methods handles evidence availability and validation", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 903)
  obj <- as_aegis(aegis_example, deconv = deconv)

  expect_error(score_methods(obj), "No method-evidence sources")

  obj <- audit_basic(obj)
  obj <- score_methods(obj)
  expect_true(is.data.frame(obj$consensus$method_evidence))
  expect_true(all(is.finite(obj$consensus$method_evidence$overall_score_raw)))
})

test_that("rank_methods supports mean_rank and rra modes", {
  obj <- make_ranking_ready_obj(seed = 904)
  obj <- score_methods(obj)

  obj_mean <- rank_methods(obj, method = "mean_rank")
  rk_mean <- obj_mean$consensus$method_ranking
  expect_true(is.data.frame(rk_mean))
  expect_true(all(c("method", "overall_rank", "overall_score", "recommendation") %in% colnames(rk_mean)))

  obj_rra <- rank_methods(obj, method = "rra")
  rk_rra <- obj_rra$consensus$method_ranking
  expect_true(is.data.frame(rk_rra))
  expect_true(all(c("method", "overall_rank", "overall_score", "recommendation") %in% colnames(rk_rra)))
  expect_true(all(rk_rra$recommendation %in% c("preferred", "acceptable", "use_with_caution")))
})

test_that("compute_consensus supports mean/weighted/trimmed_mean and top_n", {
  obj <- make_ranking_ready_obj(seed = 905)
  obj <- score_methods(obj)
  obj <- rank_methods(obj, method = "mean_rank")

  obj_mean <- compute_consensus(obj, strategy = "mean")
  expect_identical(obj_mean$consensus$result$strategy_used, "mean")
  expect_true(is.matrix(obj_mean$consensus$result$consensus_matrix))
  expect_true(is.matrix(obj_mean$consensus$result$method_disagreement))
  expect_true(is.data.frame(obj_mean$consensus$result$spot_confidence))

  obj_weighted <- compute_consensus(obj, strategy = "weighted")
  expect_identical(obj_weighted$consensus$result$strategy_used, "weighted")
  expect_true("methods_used" %in% names(obj_weighted$consensus$result))
  expect_true("method_weights" %in% names(obj_weighted$consensus$result))
  expect_equal(sum(obj_weighted$consensus$result$method_weights$weight), 1, tolerance = 1e-8)

  obj_trimmed <- compute_consensus(obj, strategy = "trimmed_mean")
  expect_true(obj_trimmed$consensus$result$strategy_used %in% c("trimmed_mean", "mean"))

  obj_top2 <- compute_consensus(obj, strategy = "weighted", top_n = 2)
  expect_equal(length(obj_top2$consensus$result$methods_used), 2)
  expect_equal(obj_top2$consensus$result$top_n, 2)
})
