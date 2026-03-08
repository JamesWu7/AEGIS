test_that("as_aegis constructs a valid object", {
  skip_if_not(exists("aegis_example", where = asNamespace("AEGIS")))
  data("aegis_example", package = "AEGIS")

  deconv <- simulate_deconv_results(aegis_example, seed = 101)
  obj <- as_aegis(aegis_example, deconv, markers = aegis_default_markers())

  expect_s3_class(obj, "aegis")
  expect_true(is.list(obj$deconv))
  expect_setequal(names(obj$deconv), c("RCTD", "SPOTlight", "cell2location"))
})
