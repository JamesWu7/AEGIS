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

  spots <- colnames(seu)
  if (is.null(spots) || anyDuplicated(spots)) {
    stop("Seurat spot names (`colnames(seu)`) must be present and unique.", call. = FALSE)
  }
  if (length(spots) == 0L) {
    stop("Seurat object must contain at least one spot.", call. = FALSE)
  }

  deconv_aligned <- validate_deconv_list(deconv, seurat_spots = spots)

  if (!is.null(markers)) {
    if (!is.list(markers)) {
      stop("`markers` must be NULL or a named list.", call. = FALSE)
    }
    marker_names <- names(markers)
    if (length(markers) > 0L && (is.null(marker_names) || any(trimws(marker_names) == ""))) {
      stop("`markers`, when supplied, must be a named list.", call. = FALSE)
    }
  }

  if (is.null(meta)) {
    meta <- list()
  } else if (is.list(meta)) {
    meta <- meta
  } else if (is.atomic(meta) || is.data.frame(meta)) {
    meta <- as.list(meta)
  } else {
    stop("`meta` must be NULL, a list, or coercible to a list.", call. = FALSE)
  }

  out <- list(
    seu = seu,
    deconv = deconv_aligned,
    markers = markers,
    audit = list(),
    consensus = list(),
    meta = meta
  )

  class(out) <- "aegis"
  out
}
