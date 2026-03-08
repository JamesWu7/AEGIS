# Internal utility helpers

#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
assert_is_seurat <- function(seu, arg = "seu") {
  if (!inherits(seu, "Seurat")) {
    stop(sprintf("`%s` must be a Seurat object.", arg), call. = FALSE)
  }
}

#' @keywords internal
assert_is_aegis <- function(x, arg = "x") {
  if (!inherits(x, "aegis")) {
    stop(sprintf("`%s` must be an aegis object. Use as_aegis() first.", arg), call. = FALSE)
  }
}

#' @keywords internal
coerce_to_numeric_matrix <- function(obj, method_name) {
  if (is.data.frame(obj)) {
    obj <- as.matrix(obj)
  }

  if (!is.matrix(obj)) {
    stop(
      sprintf("Method '%s': value must be coercible to matrix/data.frame.", method_name),
      call. = FALSE
    )
  }

  if (is.null(rownames(obj))) {
    stop(sprintf("Method '%s': row names are required.", method_name), call. = FALSE)
  }

  if (anyNA(rownames(obj)) || any(trimws(rownames(obj)) == "")) {
    stop(sprintf("Method '%s': missing/empty row names are not allowed.", method_name), call. = FALSE)
  }

  if (anyDuplicated(rownames(obj))) {
    stop(sprintf("Method '%s': duplicated row names are not allowed.", method_name), call. = FALSE)
  }

  suppressWarnings(storage.mode(obj) <- "double")
  if (!is.numeric(obj)) {
    stop(sprintf("Method '%s': matrix must be numeric.", method_name), call. = FALSE)
  }
  if (all(is.na(obj))) {
    stop(sprintf("Method '%s': matrix coercion resulted in all NA values; numeric input is required.", method_name), call. = FALSE)
  }
  if (anyNA(obj)) {
    stop(sprintf("Method '%s': matrix contains NA values; numeric, non-missing values are required.", method_name), call. = FALSE)
  }

  obj
}

#' @keywords internal
align_matrix_to_spots <- function(mat, spots, method_name) {
  missing_spots <- setdiff(spots, rownames(mat))
  if (length(missing_spots) > 0) {
    preview <- paste(utils::head(missing_spots, 5L), collapse = ", ")
    stop(
      sprintf(
        "Method '%s': %d Seurat spots are missing from row names (e.g. %s).",
        method_name,
        length(missing_spots),
        preview
      ),
      call. = FALSE
    )
  }

  extra_spots <- setdiff(rownames(mat), spots)
  if (length(extra_spots) > 0) {
    warning(
      sprintf(
        "Method '%s': dropping %d rows not present in Seurat object.",
        method_name,
        length(extra_spots)
      ),
      call. = FALSE
    )
    mat <- mat[setdiff(rownames(mat), extra_spots), , drop = FALSE]
  }

  if (!identical(rownames(mat), spots)) {
    warning(
      sprintf("Method '%s': reordering rows to match Seurat spot order.", method_name),
      call. = FALSE
    )
    mat <- mat[spots, , drop = FALSE]
  }

  mat
}

#' @keywords internal
get_spatial_coordinates <- function(seu) {
  assert_is_seurat(seu)

  image_names <- names(seu@images)
  if (length(image_names) == 0) {
    stop("No spatial image found in Seurat object (`seu@images` is empty).", call. = FALSE)
  }

  coords <- Seurat::GetTissueCoordinates(seu, image = image_names[[1L]])

  cn <- colnames(coords)
  if (!all(c("imagerow", "imagecol") %in% cn)) {
    if (all(c("x", "y") %in% cn)) {
      coords$imagecol <- coords$x
      coords$imagerow <- coords$y
    } else if (all(c("row", "col") %in% cn)) {
      coords$imagecol <- coords$col
      coords$imagerow <- coords$row
    }
  }

  if (!all(c("imagerow", "imagecol") %in% colnames(coords))) {
    stop("Spatial coordinates must include `imagerow` and `imagecol` (or x/y).", call. = FALSE)
  }

  if (!is.null(rownames(coords)) && all(nzchar(rownames(coords)))) {
    coords$spot <- rownames(coords)
  } else if ("cell" %in% colnames(coords)) {
    coords$spot <- as.character(coords$cell)
  } else {
    stop("Could not determine spot identifiers from tissue coordinates.", call. = FALSE)
  }
  coords
}

#' @keywords internal
method_list <- function(x) {
  names(x$deconv)
}

#' @keywords internal
shared_celltypes_all_methods <- function(deconv_list) {
  Reduce(intersect, lapply(deconv_list, colnames))
}

#' @keywords internal
row_entropy <- function(mat, eps = 1e-12) {
  p <- pmax(mat, eps)
  -rowSums(p * log(p))
}

#' @keywords internal
nearest_neighbors <- function(coords, k = 6L) {
  n <- nrow(coords)
  if (n <= 1L) {
    return(vector("list", n))
  }

  k <- max(1L, min(as.integer(k), n - 1L))
  dmat <- as.matrix(stats::dist(coords[, c("imagecol", "imagerow")]))
  diag(dmat) <- Inf

  lapply(seq_len(n), function(i) {
    order(dmat[i, ], decreasing = FALSE)[seq_len(k)]
  })
}

#' @keywords internal
resolve_template <- function() {
  pkg_template <- system.file("templates", "aegis_report.Rmd", package = "AEGIS")
  if (nzchar(pkg_template) && file.exists(pkg_template)) {
    return(pkg_template)
  }

  local_template <- file.path("inst", "templates", "aegis_report.Rmd")
  if (file.exists(local_template)) {
    return(normalizePath(local_template, winslash = "/", mustWork = TRUE))
  }

  stop("Could not locate report template `inst/templates/aegis_report.Rmd`.", call. = FALSE)
}

#' @keywords internal
resolve_path_in_dir <- function(data_dir, filename, required = TRUE) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path) && required) {
    stop(sprintf("Required file not found: %s", normalizePath(path, winslash = "/", mustWork = FALSE)), call. = FALSE)
  }
  path
}

utils::globalVariables(c(".data"))
