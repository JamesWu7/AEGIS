make_p9_seu <- function(spots = c("AAAC-1", "AAAC-2", "AAAC-3", "AAAC-4")) {
  counts <- Matrix::Matrix(
    matrix(c(2, 1, 3, 4, 5, 2, 1, 3), nrow = 2, ncol = length(spots), byrow = TRUE),
    sparse = TRUE
  )
  rownames(counts) <- c("g1", "g2")
  colnames(counts) <- spots
  Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")
}

make_runner_matrix <- function(spots, shift = 0) {
  n <- length(spots)
  mat <- cbind(
    B_cell = rep(c(0.50, 0.20, 0.10, 0.30), length.out = n) + shift,
    T_cell = rep(c(0.30, 0.50, 0.40, 0.30), length.out = n),
    Myeloid = rep(c(0.20, 0.30, 0.50, 0.40), length.out = n) - shift
  )
  rownames(mat) <- spots
  rs <- rowSums(mat)
  mat <- mat / rs
  mat
}

make_override_runner <- function(shift = 0) {
  function(seu, reference, normalize = TRUE, strict = TRUE, ...) {
    mat <- make_runner_matrix(colnames(seu), shift = shift)
    if (!isTRUE(normalize)) {
      mat <- mat * 100
    }
    mat
  }
}

test_that("get_supported_methods exposes registry metadata", {
  reg <- get_supported_methods()
  expect_true(is.data.frame(reg))
  expect_true(all(c(
    "method_name", "support_mode", "can_run_in_r", "can_run_in_python",
    "requires_reference", "adapter_reader", "reader_function", "runner_function", "dependency_type"
  ) %in% colnames(reg)))
  expect_true(all(c(
    "RCTD", "SPOTlight", "cell2location", "CARD", "SpatialDWLS",
    "stereoscope", "DestVI", "Tangram", "STdeconvolve", "DSTG", "STRIDE"
  ) %in% reg$method_name))
  expect_true(all(reg$support_mode %in% c("run_and_import_r", "run_and_import_python", "import_only", "experimental")))

  ns <- getNamespaceExports("AEGIS")
  reader_fns <- unique(stats::na.omit(reg$reader_function))
  runner_fns <- unique(stats::na.omit(reg$runner_function))
  expect_true(all(reader_fns %in% ns))
  expect_true(all(runner_fns %in% ns))
})

test_that("run_deconvolution executes selected R methods with stable output contract", {
  seu <- make_p9_seu()
  ref <- matrix(c(1, 2, 3, 4), nrow = 2)

  res <- run_deconvolution(
    seu = seu,
    reference = ref,
    methods = c("spotlight", "CARD"),
    runner_overrides = list(
      SPOTlight = make_override_runner(shift = 0.00),
      CARD = make_override_runner(shift = 0.02)
    )
  )

  expect_true(is.list(res))
  expect_true(all(c("seu", "deconv", "methods_run", "methods_skipped", "messages") %in% names(res)))
  expect_identical(sort(res$methods_run), sort(c("SPOTlight", "CARD")))
  expect_true(is.matrix(res$deconv$SPOTlight))
  expect_true(is.matrix(res$deconv$CARD))
  expect_identical(rownames(res$deconv$SPOTlight), colnames(seu))
})

test_that("run_deconvolution handles unsupported/import-only/python-disabled requests clearly", {
  seu <- make_p9_seu()
  ref <- matrix(c(1, 2, 3, 4), nrow = 2)

  expect_error(
    run_deconvolution(
      seu = seu,
      reference = ref,
      methods = "NotAMethod",
      runner_overrides = list(SPOTlight = make_override_runner())
    ),
    "Unsupported method"
  )

  expect_error(
    run_deconvolution(
      seu = seu,
      reference = ref,
      methods = "STdeconvolve",
      strict = TRUE
    ),
    "import_only"
  )

  res_import_only <- run_deconvolution(
    seu = seu,
    reference = ref,
    methods = c("STdeconvolve", "SPOTlight"),
    strict = FALSE,
    runner_overrides = list(SPOTlight = make_override_runner())
  )
  expect_true("SPOTlight" %in% res_import_only$methods_run)
  expect_true("STdeconvolve" %in% res_import_only$methods_skipped$method)

  expect_error(
    run_deconvolution(
      seu = seu,
      reference = ref,
      methods = "cell2location",
      use_python = FALSE,
      strict = TRUE
    ),
    "use_python = FALSE"
  )

  res_py_skip <- run_deconvolution(
    seu = seu,
    reference = ref,
    methods = c("cell2location", "CARD"),
    use_python = FALSE,
    strict = FALSE,
    runner_overrides = list(CARD = make_override_runner())
  )
  expect_true("CARD" %in% res_py_skip$methods_run)
  expect_true("cell2location" %in% res_py_skip$methods_skipped$method)
})

test_that("method-specific runners validate inputs and standardize outputs", {
  seu <- make_p9_seu()
  ref <- matrix(c(1, 2, 3, 4), nrow = 2)
  spots <- colnames(seu)

  result_df <- data.frame(
    spot = spots,
    B_cell = c(0.5, 0.2, 0.1, 0.3),
    T_cell = c(0.3, 0.5, 0.4, 0.3),
    Myeloid = c(0.2, 0.3, 0.5, 0.4),
    check.names = FALSE
  )

  m_sp <- run_spotlight(seu = seu, reference = ref, result = result_df)
  expect_true(is.matrix(m_sp))
  expect_identical(rownames(m_sp), spots)
  expect_equal(unname(rowSums(m_sp)), rep(1, length(spots)), tolerance = 1e-8)

  m_rctd <- run_rctd(seu = seu, reference = ref, result = list(weights = result_df))
  expect_true(is.matrix(m_rctd))
  expect_identical(rownames(m_rctd), spots)

  expect_error(
    run_card(seu = seu, reference = NULL, result = result_df),
    "requires a non-NULL `reference`"
  )
})

test_that("run_aegis_full performs one-shot deconvolution handoff", {
  data("aegis_example", package = "AEGIS")
  seu <- aegis_example
  ref <- matrix(c(1, 2, 3, 4), nrow = 2)

  obj <- run_aegis_full(
    seu = seu,
    reference = ref,
    methods = c("SPOTlight", "CARD"),
    markers = aegis_default_markers(),
    runner_overrides = list(
      SPOTlight = make_override_runner(shift = 0.00),
      CARD = make_override_runner(shift = 0.02)
    )
  )

  expect_s3_class(obj, "aegis")
  expect_true(!is.null(obj$audit$basic))
  expect_true(!is.null(obj$consensus$result))
  expect_true(all(c("SPOTlight", "CARD") %in% obj$meta$deconvolution_run$methods_run))
})
