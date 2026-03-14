fixture_p8 <- function(name) {
  testthat::test_path("fixtures", name)
}

make_small_seu <- function(spots = c("AAAC-1", "AAAC-2", "AAAC-3")) {
  counts <- Matrix::Matrix(
    matrix(
      c(1, 3, 2, 2, 1, 4),
      nrow = 2,
      ncol = length(spots),
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- c("g1", "g2")
  colnames(counts) <- spots
  Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")
}

test_that("P8 import adapters return valid matrices", {
  importers <- list(
    CARD = function() read_card(fixture_p8("card_like.csv")),
    SpatialDWLS = function() read_spatialdwls(fixture_p8("spatialdwls_like.tsv")),
    stereoscope = function() read_stereoscope(fixture_p8("stereoscope_like.csv")),
    DestVI = function() read_destvi(fixture_p8("destvi_like.csv")),
    Tangram = function() read_tangram(fixture_p8("tangram_like.tsv")),
    STdeconvolve = function() read_stdeconvolve(fixture_p8("stdeconvolve_like.csv")),
    DSTG = function() read_dstg(fixture_p8("dstg_like.csv")),
    STRIDE = function() read_stride(fixture_p8("stride_like.csv"))
  )

  for (nm in names(importers)) {
    m <- importers[[nm]]()
    expect_true(is.matrix(m), info = nm)
    expect_true(is.numeric(m), info = nm)
    expect_identical(rownames(m), c("AAAC-1", "AAAC-2", "AAAC-3"), info = nm)
    expect_false(anyNA(m), info = nm)
    expect_false(any(m < 0), info = nm)
    expect_equal(unname(rowSums(m)), c(1, 1, 1), tolerance = 1e-8, info = nm)
  }
})

test_that("generic importer works and method metadata is attached", {
  m <- read_deconv_table(
    path = fixture_p8("card_like.csv"),
    method = "CARD",
    normalize = TRUE
  )
  expect_true(is.matrix(m))
  meta <- attr(m, "aegis_method_metadata")
  expect_true(is.list(meta))
  expect_identical(meta$method_name, "CARD")
  expect_true(isTRUE(meta$normalized))
})

test_that("normalize = FALSE preserves abundance-style values", {
  m_raw <- read_destvi(fixture_p8("destvi_like.csv"), normalize = FALSE)
  expect_equal(unname(m_raw["AAAC-1", "B_cell"]), 20)
  expect_equal(unname(m_raw["AAAC-3", "Myeloid"]), 20)

  m_norm <- read_destvi(fixture_p8("destvi_like.csv"), normalize = TRUE)
  expect_equal(unname(rowSums(m_norm)), c(1, 1, 1), tolerance = 1e-8)
})

test_that("malformed and ambiguous fixtures fail clearly", {
  expect_error(read_card(fixture_p8("duplicated_spots.csv")), "duplicated spot identifiers")
  expect_error(read_card(fixture_p8("non_numeric_abundance.csv")), "no numeric cell-type columns|contains NA values|must be numeric")
  expect_error(read_spatialdwls(fixture_p8("metadata_only_numeric.csv")), "no numeric cell-type columns")
  expect_error(read_stdeconvolve(fixture_p8("stdeconvolve_malformed.csv")), "no numeric cell-type columns")
  expect_error(read_rctd(fixture_p8("ambiguous_transpose.csv"), strict = TRUE), "no numeric cell-type columns|could not infer spot identifiers")
  expect_error(read_stride(fixture_p8("stride_topic_only.csv"), strict = TRUE), "topic-only")
})

test_that("zero-only rows remain safe during normalization", {
  m <- read_deconv_table(fixture_p8("zero_only_rows.csv"), method = "generic", normalize = TRUE)
  expect_equal(unname(rowSums(m)), c(0, 1, 1), tolerance = 1e-8)
})

test_that("P8 imported results are directly usable in downstream AEGIS", {
  seu <- make_small_seu()
  deconv <- list(
    CARD = read_card(fixture_p8("card_like.csv")),
    SpatialDWLS = read_spatialdwls(fixture_p8("spatialdwls_like.tsv")),
    stereoscope = read_stereoscope(fixture_p8("stereoscope_like.csv")),
    DestVI = read_destvi(fixture_p8("destvi_like.csv"))
  )

  obj <- as_aegis(seu, deconv = deconv)
  obj <- audit_basic(obj)
  expect_true(inherits(obj, "aegis"))
  expect_true(!is.null(obj$audit$basic$summary))
})

test_that("representative P8 imported methods support compare/consensus", {
  seu <- make_small_seu()
  deconv <- list(
    CARD = read_card(fixture_p8("card_like.csv")),
    Tangram = read_tangram(fixture_p8("tangram_like.tsv")),
    DestVI = read_destvi(fixture_p8("destvi_like.csv"))
  )

  obj <- as_aegis(seu, deconv = deconv)
  obj <- compare_methods(obj)
  obj <- compute_consensus(obj)

  expect_true(!is.null(obj$consensus$comparison))
  expect_true(!is.null(obj$consensus$result$consensus_matrix))
})
