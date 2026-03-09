#' Align a deconvolution matrix to Seurat spot order
#'
#' @param deconv A deconvolution matrix/data.frame with spot row names.
#' @param seu A Seurat object.
#' @param method_name Method label used in error/warning messages.
#'
#' @return A numeric matrix aligned to `colnames(seu)`.
#' @export
align_deconv_to_seurat <- function(deconv, seu, method_name = "deconv") {
  assert_is_seurat(seu, "seu")
  spots <- colnames(seu)
  if (is.null(spots) || length(spots) == 0L || anyDuplicated(spots)) {
    stop("Seurat object must contain unique spot names in `colnames(seu)`.", call. = FALSE)
  }
  mat <- coerce_to_numeric_matrix(deconv, method_name = method_name)
  align_matrix_to_spots(mat, spots, method_name = method_name)
}

#' Read RCTD deconvolution outputs
#'
#' Reads exported RCTD result files and standardizes them to a spot-by-celltype
#' numeric matrix usable by `as_aegis()`.
#'
#' @param path File path to an exported RCTD result (csv/tsv/txt/rds).
#' @param type Input type: `auto`, `weights`, `proportions`, `results_df`, `rds`.
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_rctd <- function(
    path,
    type = c("auto", "weights", "proportions", "results_df", "rds"),
    spot_col = NULL,
    normalize = TRUE,
    strict = TRUE) {
  type <- match.arg(type)
  obj <- read_deconv_input(path = path, type = type, method = "RCTD", strict = strict)

  mat <- standardize_deconv_matrix(
    obj = obj,
    spot_col = spot_col,
    strict = strict,
    method = "RCTD"
  )
  if (isTRUE(normalize)) {
    mat <- normalize_deconv_rows(mat)
  }
  mat
}

#' Read SPOTlight deconvolution outputs
#'
#' Reads exported SPOTlight result files and standardizes them to a
#' spot-by-celltype numeric matrix usable by `as_aegis()`.
#'
#' @param path File path to an exported SPOTlight result (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_spotlight <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  obj <- read_deconv_input(path = path, type = "auto", method = "SPOTlight", strict = strict)
  mat <- standardize_deconv_matrix(
    obj = obj,
    spot_col = spot_col,
    strict = strict,
    method = "SPOTlight"
  )
  if (isTRUE(normalize)) {
    mat <- normalize_deconv_rows(mat)
  }
  mat
}

#' Read cell2location deconvolution outputs
#'
#' Reads exported cell2location abundance/proportion tables and standardizes them
#' to a spot-by-celltype numeric matrix usable by `as_aegis()`.
#'
#' @param path File path to an exported cell2location result (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, abundances are row-normalized to proportions.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_cell2location <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  obj <- read_deconv_input(path = path, type = "auto", method = "cell2location", strict = strict)
  mat <- standardize_deconv_matrix(
    obj = obj,
    spot_col = spot_col,
    strict = strict,
    method = "cell2location"
  )
  if (isTRUE(normalize)) {
    mat <- normalize_deconv_rows(mat)
  }
  mat
}

#' @keywords internal
validate_deconv_list <- function(deconv, seurat_spots = NULL) {
  if (!is.list(deconv) || length(deconv) == 0L) {
    stop("`deconv` must be a non-empty named list.", call. = FALSE)
  }

  method_names <- names(deconv)
  if (is.null(method_names)) {
    stop("`deconv` must be a named list with non-empty method names.", call. = FALSE)
  }
  if (anyNA(method_names) || any(trimws(method_names) == "")) {
    stop("`deconv` contains empty method names. Every method must be named.", call. = FALSE)
  }
  if (anyDuplicated(method_names)) {
    stop("`deconv` method names must be unique.", call. = FALSE)
  }

  out <- lapply(seq_along(deconv), function(i) {
    method <- method_names[[i]]
    mat <- coerce_to_numeric_matrix(deconv[[i]], method_name = method)
    if (!is.null(seurat_spots)) {
      mat <- align_matrix_to_spots(mat, spots = seurat_spots, method_name = method)
    }
    mat
  })
  names(out) <- method_names
  out
}

#' @keywords internal
read_deconv_input <- function(path, type = c("auto", "rds", "weights", "proportions", "results_df"), method = "deconv", strict = TRUE) {
  type <- match.arg(type)
  if (!is.character(path) || length(path) != 1L || is.na(path) || trimws(path) == "") {
    stop("`path` must be a single non-empty file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("%s input file does not exist: %s", method, normalizePath(path, winslash = "/", mustWork = FALSE)), call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  use_rds <- identical(type, "rds") || (identical(type, "auto") && ext == "rds")

  if (use_rds) {
    obj <- readRDS(path)
    return(extract_matrix_like_from_object(obj, type = type, method = method, strict = strict))
  }

  if (ext == "rds") {
    obj <- readRDS(path)
    return(extract_matrix_like_from_object(obj, type = type, method = method, strict = strict))
  }

  read_delimited_deconv(path)
}

#' @keywords internal
read_delimited_deconv <- function(path) {
  first_line <- readLines(path, n = 1L, warn = FALSE)
  sep <- ","
  if (length(first_line) == 1L) {
    n_tab <- lengths(regmatches(first_line, gregexpr("\t", first_line, fixed = TRUE)))
    n_com <- lengths(regmatches(first_line, gregexpr(",", first_line, fixed = TRUE)))
    if (n_tab > n_com) {
      sep <- "\t"
    } else if (n_tab == 0L && n_com == 0L) {
      sep <- ""
    }
  }

  out <- utils::read.table(
    path,
    header = TRUE,
    sep = sep,
    quote = "\"",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    fill = TRUE
  )
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' @keywords internal
extract_matrix_like_from_object <- function(obj, type = "auto", method = "deconv", strict = TRUE) {
  if (is.matrix(obj) || is.data.frame(obj)) {
    return(obj)
  }
  if (!is.list(obj)) {
    stop(sprintf("%s RDS must contain a matrix/data.frame or a list with one.", method), call. = FALSE)
  }

  preferred <- switch(type,
    weights = c("weights", "norm_weights"),
    proportions = c("proportions", "weights", "norm_weights"),
    results_df = c("results_df", "results", "weights"),
    auto = c("weights", "proportions", "results_df", "results", "norm_weights"),
    rds = c("weights", "proportions", "results_df", "results", "norm_weights"),
    c("weights", "proportions", "results_df", "results", "norm_weights")
  )

  for (nm in preferred) {
    if (nm %in% names(obj) && (is.matrix(obj[[nm]]) || is.data.frame(obj[[nm]]))) {
      return(obj[[nm]])
    }
  }

  matrix_candidates <- names(obj)[vapply(obj, function(x) is.matrix(x) || is.data.frame(x), logical(1))]
  if (length(matrix_candidates) == 1L) {
    return(obj[[matrix_candidates]])
  }
  if (length(matrix_candidates) > 1L && isTRUE(strict)) {
    stop(
      sprintf(
        "%s RDS contains multiple matrix-like entries (%s). Please set `type` or provide a cleaner exported table.",
        method,
        paste(matrix_candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(matrix_candidates) > 1L) {
    return(obj[[matrix_candidates[[1L]]]])
  }

  stop(sprintf("Could not extract a matrix/data.frame from %s RDS input.", method), call. = FALSE)
}

#' @keywords internal
standardize_deconv_matrix <- function(obj, spot_col = NULL, strict = TRUE, method = "deconv") {
  if (is.data.frame(obj)) {
    df <- as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.matrix(obj)) {
    df <- as.data.frame(obj, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop(sprintf("%s input must be matrix/data.frame after reading.", method), call. = FALSE)
  }

  if (nrow(df) == 0L || ncol(df) == 0L) {
    stop(sprintf("%s table is empty.", method), call. = FALSE)
  }

  spots <- NULL
  if (!is.null(spot_col)) {
    if (!(spot_col %in% colnames(df))) {
      stop(sprintf("%s: `spot_col`='%s' not found in input columns.", method, spot_col), call. = FALSE)
    }
    spots <- as.character(df[[spot_col]])
    df[[spot_col]] <- NULL
  } else if (!is.null(rownames(df)) && !identical(rownames(df), as.character(seq_len(nrow(df))))) {
    spots <- rownames(df)
  } else {
    guessed_spot <- guess_spot_column(df)
    if (!is.null(guessed_spot)) {
      spots <- as.character(df[[guessed_spot]])
      df[[guessed_spot]] <- NULL
    }
  }

  if (is.null(spots)) {
    stop(
      sprintf("%s: could not infer spot identifiers. Provide `spot_col` or include row names / a spot barcode column.", method),
      call. = FALSE
    )
  }
  if (anyNA(spots) || any(trimws(spots) == "")) {
    stop(sprintf("%s: spot identifiers contain missing or empty values.", method), call. = FALSE)
  }
  if (anyDuplicated(spots)) {
    stop(sprintf("%s: duplicated spot identifiers detected.", method), call. = FALSE)
  }

  cell_cols <- guess_celltype_columns(df, strict = strict)
  if (length(cell_cols) == 0L) {
    stop(sprintf("%s: no numeric cell-type columns could be identified.", method), call. = FALSE)
  }

  dropped <- setdiff(colnames(df), cell_cols)
  if (length(dropped) > 0L && isTRUE(strict)) {
    numeric_dropped <- dropped[vapply(df[dropped], is.numeric, logical(1))]
    if (length(numeric_dropped) > 0L) {
      warning(
        sprintf("%s: dropped likely metadata numeric columns: %s", method, paste(numeric_dropped, collapse = ", ")),
        call. = FALSE
      )
    }
  }

  mat <- as.matrix(df[, cell_cols, drop = FALSE])
  rownames(mat) <- spots
  colnames(mat) <- cell_cols

  mat <- maybe_transpose_deconv_matrix(mat, strict = strict, method = method)
  mat <- assert_numeric_deconv(mat, method = method)

  if (any(mat < 0)) {
    if (isTRUE(strict)) {
      stop(sprintf("%s matrix contains negative values; please provide non-negative abundances/proportions.", method), call. = FALSE)
    }
    warning(sprintf("%s matrix contains negative values.", method), call. = FALSE)
  }

  mat
}

#' @keywords internal
guess_spot_column <- function(df) {
  stopifnot(is.data.frame(df))
  if (ncol(df) == 0L) return(NULL)

  norm_names <- tolower(gsub("[^a-z0-9]+", "", colnames(df)))
  candidates <- c("spot", "spots", "spotid", "barcode", "barcodes", "cell", "cellid", "id", "rownames", "rowname")
  hit <- which(norm_names %in% candidates)
  if (length(hit) >= 1L) {
    cand <- colnames(df)[hit[[1L]]]
    vals <- as.character(df[[cand]])
    numeric_like <- suppressWarnings(!any(is.na(as.numeric(vals))))
    if (!numeric_like) {
      return(cand)
    }
  }

  first_name <- colnames(df)[[1L]]
  first_norm <- norm_names[[1L]]
  first_vals <- as.character(df[[1L]])
  unique_ok <- !anyDuplicated(first_vals) && all(nzchar(first_vals))
  numeric_like <- suppressWarnings(!any(is.na(as.numeric(first_vals))))

  if (unique_ok && !numeric_like && (first_norm %in% c("x", "unnamed0", "v1", "...1") || !numeric_like)) {
    return(first_name)
  }

  NULL
}

#' @keywords internal
guess_celltype_columns <- function(df, strict = TRUE) {
  stopifnot(is.data.frame(df))
  if (ncol(df) == 0L) return(character(0))

  meta_like <- c(
    "sample", "sampleid", "sample_id", "x", "y", "row", "col",
    "imagerow", "imagecol", "pxlrowinfullres", "pxlcolinfullres",
    "confidence", "cluster", "region", "slice", "tissue", "library", "batch"
  )
  norm_names <- tolower(gsub("[^a-z0-9]+", "", colnames(df)))

  is_numeric_like <- vapply(df, function(v) {
    if (is.numeric(v)) return(TRUE)
    vv <- suppressWarnings(as.numeric(v))
    sum(!is.na(vv)) >= max(2L, floor(length(v) * 0.8))
  }, logical(1))

  is_meta_name <- norm_names %in% meta_like
  keep <- is_numeric_like & !is_meta_name

  if (sum(keep) == 0L && !isTRUE(strict)) {
    keep <- is_numeric_like
  }
  colnames(df)[keep]
}

#' @keywords internal
maybe_transpose_deconv_matrix <- function(mat, strict = TRUE, method = "deconv", seurat_spots = NULL) {
  if (!is.matrix(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
    stop(sprintf("%s matrix is empty or invalid.", method), call. = FALSE)
  }
  rn <- rownames(mat)
  cn <- colnames(mat)
  if (is.null(cn)) {
    stop(sprintf("%s matrix must include column names.", method), call. = FALSE)
  }

  is_spot_like <- function(x) {
    if (is.null(x)) return(rep(FALSE, 0L))
    x <- as.character(x)
    grepl("-[0-9]+$", x) |
      grepl("_[0-9]+$", x) |
      grepl("^[ACGTN]+-[0-9]+$", x, ignore.case = TRUE) |
      (nchar(x) >= 14L & grepl("[0-9]", x))
  }
  is_celltype_like <- function(x) {
    if (is.null(x)) return(rep(FALSE, 0L))
    x <- as.character(x)
    has_letter <- grepl("[A-Za-z]", x)
    not_pure_numeric <- !grepl("^[0-9]+$", x)
    not_barcode <- !grepl("-[0-9]+$", x) & !grepl("_[0-9]+$", x)
    has_letter & not_pure_numeric & not_barcode
  }

  if (!is.null(seurat_spots)) {
    row_hit <- if (is.null(rn)) 0 else mean(rn %in% seurat_spots)
    col_hit <- if (is.null(cn)) 0 else mean(cn %in% seurat_spots)
    if (col_hit > row_hit && col_hit >= 0.6 && row_hit < 0.2) {
      mat <- t(mat)
      return(mat)
    }
    if (row_hit > col_hit && row_hit >= 0.6) {
      return(mat)
    }
    if (isTRUE(strict) && row_hit > 0.3 && col_hit > 0.3) {
      stop(sprintf("%s matrix orientation is ambiguous against Seurat spot names.", method), call. = FALSE)
    }
    return(mat)
  }

  row_spot <- if (is.null(rn)) 0 else mean(is_spot_like(rn))
  col_spot <- mean(is_spot_like(cn))
  row_ct <- if (is.null(rn)) 0 else mean(is_celltype_like(rn))
  col_ct <- mean(is_celltype_like(cn))

  should_transpose <- FALSE
  if (!is.null(rn) &&
    nrow(mat) <= 80L &&
    ncol(mat) >= nrow(mat) &&
    row_ct >= 0.7 &&
    col_spot >= 0.5 &&
    row_spot <= 0.4) {
    should_transpose <- TRUE
  }

  if (isTRUE(should_transpose)) {
    mat <- t(mat)
  } else if (isTRUE(strict) && !is.null(rn) && row_spot >= 0.6 && col_spot >= 0.6) {
    stop(sprintf("%s matrix orientation is ambiguous; please provide clearer row/column labels.", method), call. = FALSE)
  }

  mat
}

#' @keywords internal
normalize_deconv_rows <- function(mat) {
  mat <- assert_numeric_deconv(mat, method = "deconv")
  rs <- rowSums(mat)
  positive <- is.finite(rs) & rs > 0
  if (any(positive)) {
    mat[positive, ] <- mat[positive, , drop = FALSE] / rs[positive]
  }
  mat
}

#' @keywords internal
assert_numeric_deconv <- function(mat, method = "deconv") {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  suppressWarnings(storage.mode(mat) <- "double")
  if (!is.numeric(mat)) {
    stop(sprintf("%s matrix must be numeric.", method), call. = FALSE)
  }
  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    stop(sprintf("%s matrix has zero rows or columns.", method), call. = FALSE)
  }
  if (is.null(rownames(mat)) || anyNA(rownames(mat)) || any(trimws(rownames(mat)) == "")) {
    stop(sprintf("%s matrix must have non-empty rownames (spots).", method), call. = FALSE)
  }
  if (anyDuplicated(rownames(mat))) {
    stop(sprintf("%s matrix has duplicated spot rownames.", method), call. = FALSE)
  }
  if (is.null(colnames(mat)) || anyNA(colnames(mat)) || any(trimws(colnames(mat)) == "")) {
    stop(sprintf("%s matrix must have non-empty column names (cell types).", method), call. = FALSE)
  }
  if (anyDuplicated(colnames(mat))) {
    stop(sprintf("%s matrix has duplicated cell-type column names.", method), call. = FALSE)
  }
  if (anyNA(mat)) {
    stop(sprintf("%s matrix contains NA values.", method), call. = FALSE)
  }
  if (any(!is.finite(mat))) {
    stop(sprintf("%s matrix contains non-finite values.", method), call. = FALSE)
  }
  mat
}
