fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

test_that("read_rctd imports delimited and RDS inputs", {
  m_csv <- read_rctd(fixture_path("rctd_like.csv"), normalize = TRUE)
  expect_true(is.matrix(m_csv))
  expect_identical(rownames(m_csv), c("AAAC-1", "AAAC-2", "AAAC-3"))
  expect_identical(colnames(m_csv), c("B_cell", "T_cell", "Myeloid"))
  expect_equal(unname(rowSums(m_csv)[1:2]), c(1, 1), tolerance = 1e-8)
  expect_equal(unname(rowSums(m_csv)[3]), 0)

  m_rds <- read_rctd(fixture_path("rctd_like_matrix.rds"), type = "rds")
  expect_true(is.matrix(m_rds))
  expect_identical(rownames(m_rds), c("AAAC-1", "AAAC-2", "AAAC-3"))

  m_list <- read_rctd(fixture_path("rctd_like_list.rds"), type = "auto")
  expect_true(is.matrix(m_list))
  expect_identical(colnames(m_list), c("B_cell", "T_cell", "Myeloid"))
})

test_that("read_rctd handles malformed and strict ambiguity paths", {
  expect_error(
    read_rctd(fixture_path("malformed_deconv.csv")),
    "no numeric cell-type columns"
  )
  expect_error(
    read_rctd(fixture_path("no_spot_numeric.csv"), strict = TRUE),
    "could not infer spot identifiers"
  )
})

test_that("read_rctd can recover transposed-like table", {
  m <- read_rctd(fixture_path("transposed_like.csv"), normalize = FALSE)
  expect_true(is.matrix(m))
  expect_identical(rownames(m), c("AAAC-1", "AAAC-2", "AAAC-3"))
  expect_identical(colnames(m), c("B_cell", "T_cell", "Myeloid"))
})

test_that("read_spotlight imports and excludes metadata columns", {
  m <- read_spotlight(fixture_path("spotlight_like.tsv"), normalize = TRUE)
  expect_true(is.matrix(m))
  expect_identical(rownames(m), c("AAAC-1", "AAAC-2", "AAAC-3"))
  expect_identical(colnames(m), c("B_cell", "T_cell", "NK_cell"))
  expect_equal(unname(rowSums(m)), c(1, 1, 1), tolerance = 1e-8)

  expect_error(read_spotlight(fixture_path("malformed_deconv.csv")), "no numeric cell-type columns")
})

test_that("read_cell2location supports normalized and raw abundance modes", {
  m_norm <- read_cell2location(fixture_path("cell2location_like.csv"), normalize = TRUE)
  expect_true(is.matrix(m_norm))
  expect_identical(colnames(m_norm), c("B_cell", "T_cell", "Myeloid"))
  expect_equal(unname(rowSums(m_norm)[1:2]), c(1, 1), tolerance = 1e-8)
  expect_equal(unname(rowSums(m_norm)[3]), 0)

  m_raw <- read_cell2location(fixture_path("cell2location_like.csv"), normalize = FALSE)
  expect_equal(unname(m_raw["AAAC-1", "B_cell"]), 20)
  expect_error(read_cell2location(fixture_path("malformed_deconv.csv")), "no numeric cell-type columns")
})

test_that("helper functions validate and align matrices safely", {
  good <- list(
    RCTD = read_rctd(fixture_path("rctd_like.csv")),
    SPOTlight = read_spotlight(fixture_path("spotlight_like.tsv"))
  )
  out <- AEGIS:::validate_deconv_list(good)
  expect_true(is.list(out))
  expect_identical(names(out), c("RCTD", "SPOTlight"))

  bad_unnamed <- good
  names(bad_unnamed) <- NULL
  expect_error(AEGIS:::validate_deconv_list(bad_unnamed), "named list")

  counts <- Matrix::Matrix(
    matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3, byrow = TRUE),
    sparse = TRUE
  )
  rownames(counts) <- c("g1", "g2")
  colnames(counts) <- c("AAAC-1", "AAAC-2", "AAAC-3")
  seu_small <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")

  mat <- read_rctd(fixture_path("rctd_like.csv"), normalize = FALSE)
  perm <- mat[c("AAAC-2", "AAAC-1", "AAAC-3"), , drop = FALSE]
  expect_warning(
    aligned <- align_deconv_to_seurat(perm, seu_small, method_name = "RCTD"),
    "reordering rows"
  )
  expect_identical(rownames(aligned), colnames(seu_small))
})

test_that("standardization and normalization helpers behave on edge cases", {
  raw <- utils::read.csv(fixture_path("rctd_like.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  std <- AEGIS:::standardize_deconv_matrix(raw, method = "test")
  expect_identical(rownames(std), c("AAAC-1", "AAAC-2", "AAAC-3"))
  expect_identical(colnames(std), c("B_cell", "T_cell", "Myeloid"))

  m <- matrix(c(2, 2, 0, 0), nrow = 2, byrow = TRUE)
  rownames(m) <- c("s1", "s2")
  colnames(m) <- c("ct1", "ct2")
  m2 <- AEGIS:::normalize_deconv_rows(m)
  expect_equal(unname(rowSums(m2)), c(1, 0))

  keep <- matrix(c(0.6, 0.4, 0.2, 0.8), nrow = 2, byrow = TRUE)
  rownames(keep) <- c("AAAC-1", "AAAC-2")
  colnames(keep) <- c("B_cell", "T_cell")
  expect_identical(
    AEGIS:::maybe_transpose_deconv_matrix(keep, strict = TRUE, method = "test"),
    keep
  )
})

test_that("real-input workflow patterns run with imported tables", {
  data("aegis_example", package = "AEGIS")
  seu <- suppressWarnings(aegis_example[, colnames(aegis_example)[1:3]])
  spots <- colnames(seu)

  tmp_rctd <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      barcode = spots,
      B_cell = c(0.4, 0.2, 0.1),
      T_cell = c(0.4, 0.5, 0.3),
      Myeloid = c(0.2, 0.3, 0.6),
      check.names = FALSE
    ),
    tmp_rctd,
    row.names = FALSE
  )

  rctd <- read_rctd(tmp_rctd)
  obj1 <- as_aegis(seu, deconv = list(RCTD = rctd))
  obj1 <- audit_basic(obj1)
  expect_true(inherits(obj1, "aegis"))
  expect_true(!is.null(obj1$audit$basic))

  tmp_sp <- tempfile(fileext = ".tsv")
  utils::write.table(
    data.frame(
      spot_id = spots,
      B_cell = c(0.3, 0.4, 0.1),
      T_cell = c(0.5, 0.4, 0.5),
      Myeloid = c(0.2, 0.2, 0.4),
      sample = "S1",
      check.names = FALSE
    ),
    tmp_sp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  tmp_c2l <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      spot = spots,
      B_cell = c(15, 10, 2),
      T_cell = c(18, 11, 10),
      Myeloid = c(4, 3, 9),
      x = c(1, 2, 3),
      y = c(5, 6, 7),
      check.names = FALSE
    ),
    tmp_c2l,
    row.names = FALSE
  )

  spotlight <- read_spotlight(tmp_sp)
  cell2location <- read_cell2location(tmp_c2l)
  obj2 <- as_aegis(seu, deconv = list(SPOTlight = spotlight, cell2location = cell2location))
  obj2 <- audit_basic(obj2)
  obj2 <- compare_methods(obj2)
  obj2 <- compute_consensus(obj2)
  expect_true(inherits(obj2, "aegis"))
  expect_true(!is.null(obj2$consensus$comparison))
  expect_true(!is.null(obj2$consensus$result))
})
