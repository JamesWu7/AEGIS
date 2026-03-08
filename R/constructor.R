#' Construct an AEGIS object
#'
#' Creates a lightweight S3 `aegis` object that stores a Seurat object,
#' deconvolution matrices, marker definitions, and downstream analysis state.
#'
#' @param seu A Seurat object with spatial spots in columns.
#' @param deconv Named list of deconvolution matrices/data frames.
#'   Each element must be spot-by-celltype with row names aligned to Seurat spot names.
#' @param markers Optional named list of marker genes by cell type.
#' @param meta Optional list of user metadata.
#'
#' @return An object of class `aegis`.
#' @export
as_aegis <- function(seu, deconv, markers = NULL, meta = NULL) {
  assert_is_seurat(seu, "seu")

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

  spots <- colnames(seu)
  if (is.null(spots) || anyDuplicated(spots)) {
    stop("Seurat spot names (`colnames(seu)`) must be present and unique.", call. = FALSE)
  }

  deconv_aligned <- lapply(seq_along(deconv), function(i) {
    method <- method_names[[i]]
    mat <- coerce_to_numeric_matrix(deconv[[i]], method)
    align_matrix_to_spots(mat, spots, method)
  })
  names(deconv_aligned) <- method_names

  if (!is.null(markers) && !is.list(markers)) {
    stop("`markers` must be NULL or a named list.", call. = FALSE)
  }

  if (!is.null(meta) && !is.list(meta)) {
    stop("`meta` must be NULL or a list.", call. = FALSE)
  }

  out <- list(
    seu = seu,
    deconv = deconv_aligned,
    markers = markers,
    audit = list(),
    consensus = list(),
    meta = c(list(created_at = as.character(Sys.time())), meta %||% list())
  )

  class(out) <- "aegis"
  out
}
