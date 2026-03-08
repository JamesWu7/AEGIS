test_that("simulate_deconv_results returns aligned matrices", {
  data("aegis_example", package = "AEGIS")

  deconv <- simulate_deconv_results(aegis_example, seed = 7)
  expect_true(is.list(deconv))
  expect_true(length(deconv) >= 2)

  for (m in names(deconv)) {
    mat <- deconv[[m]]
    expect_true(is.matrix(mat))
    expect_identical(rownames(mat), colnames(aegis_example))
    expect_true(all(is.finite(mat)))
    expect_true(all(abs(rowSums(mat) - 1) < 1e-6))
  }
})
