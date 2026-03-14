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

#' Read a generic deconvolution result table
#'
#' Generic importer for spot-by-celltype tables exported by external methods.
#'
#' @param path File path to exported results (csv/tsv/txt/rds).
#' @param method Optional method label used in messages and metadata.
#' @param type Input type for RDS-aware parsing (`auto`, `rds`, `weights`,
#'   `proportions`, `results_df`).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_deconv_table <- function(
    path,
    method = NULL,
    type = c("auto", "rds", "weights", "proportions", "results_df"),
    spot_col = NULL,
    normalize = TRUE,
    strict = TRUE) {
  type <- match.arg(type)
  method <- if (is.null(method) || !nzchar(trimws(method))) "deconv" else as.character(method)[1L]

  obj <- read_deconv_input(path = path, type = type, method = method, strict = strict)
  mat <- standardize_deconv_matrix(obj = obj, spot_col = spot_col, strict = strict, method = method)
  if (isTRUE(normalize)) {
    mat <- normalize_deconv_rows(mat)
  }
  attach_method_metadata(
    mat,
    method_name = method,
    imported_from = tools::file_ext(path),
    normalized = isTRUE(normalize),
    original_file = normalizePath(path, winslash = "/", mustWork = FALSE)
  )
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
  read_deconv_table(path = path, method = "RCTD", type = match.arg(type), spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read SPOTlight deconvolution outputs
#'
#' @param path File path to an exported SPOTlight result (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_spotlight <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  read_deconv_table(path = path, method = "SPOTlight", spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read cell2location deconvolution outputs
#'
#' @param path File path to an exported cell2location result (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, abundances are row-normalized to proportions.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_cell2location <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  read_deconv_table(path = path, method = "cell2location", spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read CARD deconvolution outputs
#'
#' @param path File path to exported CARD results (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_card <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  read_deconv_table(path = path, method = "CARD", spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read SpatialDWLS deconvolution outputs
#'
#' @param path File path to exported SpatialDWLS results (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_spatialdwls <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  read_deconv_table(path = path, method = "SpatialDWLS", spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read stereoscope deconvolution outputs
#'
#' @param path File path to exported stereoscope results (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_stereoscope <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  read_deconv_table(path = path, method = "stereoscope", spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read DestVI deconvolution outputs
#'
#' @param path File path to exported DestVI abundance/proportion results (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_destvi <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  read_deconv_table(path = path, method = "DestVI", spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read Tangram mapping-derived composition outputs
#'
#' @param path File path to exported Tangram mapping-derived composition table.
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_tangram <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  read_deconv_table(path = path, method = "Tangram", spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read STdeconvolve outputs
#'
#' Supports de novo/latent cell-type labels (for example `topic1`, `topic2`).
#'
#' @param path File path to exported STdeconvolve results (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_stdeconvolve <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  read_deconv_table(path = path, method = "STdeconvolve", spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read DSTG deconvolution outputs
#'
#' @param path File path to exported DSTG results (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous parsing/transposition raises clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_dstg <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  read_deconv_table(path = path, method = "DSTG", spot_col = spot_col, normalize = normalize, strict = strict)
}

#' Read STRIDE deconvolution outputs
#'
#' @param path File path to exported STRIDE results (csv/tsv/txt/rds).
#' @param spot_col Optional explicit spot/barcode column name.
#' @param normalize If `TRUE`, rows are normalized to sum 1 when possible.
#' @param strict If `TRUE`, ambiguous topic-only files raise clear errors.
#'
#' @return A numeric spot-by-celltype matrix.
#' @export
read_stride <- function(path, spot_col = NULL, normalize = TRUE, strict = TRUE) {
  mat <- read_deconv_table(path = path, method = "STRIDE", spot_col = spot_col, normalize = normalize, strict = strict)
  if (isTRUE(strict)) {
    nm <- colnames(mat)
    topic_like <- grepl("^topic[_-]?[0-9]+$", tolower(nm))
    if (length(topic_like) > 0L && all(topic_like)) {
      stop(
        "STRIDE import appears topic-only (no direct cell-type composition columns). Export cell-type composition table before import.",
        call. = FALSE
      )
    }
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
    drop_idx <- which(colnames(df) == spot_col)
    if (length(drop_idx) == 0L) {
      stop(sprintf("%s: failed to remove `spot_col`='%s' from table.", method, spot_col), call. = FALSE)
    }
    spots <- as.character(df[[drop_idx[[1L]]]])
    df <- df[, -drop_idx[[1L]], drop = FALSE]
  } else if (!is.null(rownames(df)) && !identical(rownames(df), as.character(seq_len(nrow(df))))) {
    spots <- rownames(df)
  } else {
    guessed_spot <- guess_spot_column(df)
    if (!is.null(guessed_spot)) {
      drop_idx <- which(colnames(df) == guessed_spot)
      if (length(drop_idx) == 0L) {
        stop(sprintf("%s: failed to remove inferred spot column '%s'.", method, guessed_spot), call. = FALSE)
      }
      spots <- as.character(df[[drop_idx[[1L]]]])
      df <- df[, -drop_idx[[1L]], drop = FALSE]
    }
  }

  if (is.null(spots)) {
    stop(
      sprintf(
        "%s: could not infer spot identifiers. Expected row names or a barcode-like column (e.g. spot/barcode/cell/id). Detected columns: %s. Provide `spot_col` explicitly.",
        method,
        paste(colnames(df), collapse = ", ")
      ),
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
    stop(
      sprintf(
        "%s: no numeric cell-type columns could be identified. Columns seen: %s. Remove metadata columns or provide a cleaner exported table.",
        method,
        paste(colnames(df), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  dropped <- setdiff(colnames(df), cell_cols)
  if (length(dropped) > 0L && isTRUE(strict)) {
    dropped_present <- intersect(dropped, colnames(df))
    numeric_dropped <- character(0)
    if (length(dropped_present) > 0L) {
      numeric_dropped <- dropped_present[vapply(df[dropped_present], is.numeric, logical(1))]
    }
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
    "confidence", "conf", "cluster", "region", "slice", "tissue", "library", "batch",
    "pvalue", "qvalue", "fdr", "pval", "adjp", "likelihood", "weight", "score"
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

#' @keywords internal
attach_method_metadata <- function(mat, method_name, imported_from, normalized, original_file) {
  attr(mat, "aegis_method_metadata") <- list(
    method_name = method_name,
    imported_from = imported_from,
    normalized = isTRUE(normalized),
    original_file = original_file,
    parser_version = "p8"
  )
  mat
}
