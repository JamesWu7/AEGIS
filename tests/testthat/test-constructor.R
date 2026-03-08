test_that("as_aegis accepts valid inputs", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 101)

  obj <- as_aegis(aegis_example, deconv, markers = NULL, meta = NULL)
  expect_s3_class(obj, "aegis")
  expect_true(is.list(obj$deconv))
  expect_true(is.null(obj$markers))
  expect_true(is.list(obj$meta))
})

test_that("as_aegis validates malformed deconv inputs", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 102)

  bad_unnamed <- unname(deconv)
  expect_error(as_aegis(aegis_example, bad_unnamed), "named list")

  bad_non_numeric <- deconv
  bad_non_numeric[[1]] <- matrix(
    "bad",
    nrow = nrow(deconv[[1]]),
    ncol = ncol(deconv[[1]]),
    dimnames = dimnames(deconv[[1]])
  )
  expect_error(as_aegis(aegis_example, bad_non_numeric), "numeric")

  bad_no_rownames <- deconv
  rownames(bad_no_rownames[[1]]) <- NULL
  expect_error(as_aegis(aegis_example, bad_no_rownames), "row names")

  bad_mismatch <- deconv
  bad_mismatch[[1]] <- bad_mismatch[[1]][-1, , drop = FALSE]
  expect_error(as_aegis(aegis_example, bad_mismatch), "missing")
})

test_that("as_aegis reorders rows explicitly when names are permuted", {
  data("aegis_example", package = "AEGIS")
  deconv <- simulate_deconv_results(aegis_example, seed = 103)

  deconv_perm <- deconv
  perm <- sample(seq_len(nrow(deconv_perm[[1]])))
  deconv_perm[[1]] <- deconv_perm[[1]][perm, , drop = FALSE]
  expect_warning(obj <- as_aegis(aegis_example, deconv_perm), "reordering rows")
  expect_identical(rownames(obj$deconv[[1]]), colnames(aegis_example))
})
