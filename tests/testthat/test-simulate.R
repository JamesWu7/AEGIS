test_that("simulate_deconv_results returns valid default outputs", {
  data("aegis_example", package = "AEGIS")

  out <- simulate_deconv_results(aegis_example, seed = 7)
  expect_true(is.list(out))
  expect_identical(names(out), c("RCTD", "SPOTlight", "cell2location"))

  n_spots <- ncol(aegis_example)
  n_types <- length(colnames(out[[1]]))

  for (m in names(out)) {
    mat <- out[[m]]
    expect_true(is.matrix(mat))
    expect_true(is.numeric(mat))
    expect_identical(rownames(mat), colnames(aegis_example))
    expect_equal(nrow(mat), n_spots)
    expect_equal(ncol(mat), n_types)
    expect_true(all(is.finite(mat)))
    expect_true(all(mat >= 0))
    expect_true(all(abs(rowSums(mat) - 1) < 1e-6))
  }
})

test_that("simulate_deconv_results is reproducible with same seed", {
  data("aegis_example", package = "AEGIS")

  out1 <- simulate_deconv_results(aegis_example, seed = 123)
  out2 <- simulate_deconv_results(aegis_example, seed = 123)
  out3 <- simulate_deconv_results(aegis_example, seed = 124)

  expect_equal(out1, out2)
  expect_false(isTRUE(all.equal(out1[[1]], out3[[1]])))
})
