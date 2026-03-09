#' Run the standard AEGIS workflow with one function
#'
#' Convenience entry point that constructs an AEGIS object (if needed) and runs
#' the core audit/comparison pipeline. Supports:
#' \itemize{
#'   \item a Seurat object + deconvolution list (single-sample)
#'   \item a named list of Seurat objects + nested deconvolution list (multi-sample)
#'   \item an existing `aegis` or `aegis_multi` object
#' }
#'
#' @param x A Seurat object, named list of Seurat objects, `aegis`, or `aegis_multi`.
#' @param deconv Optional deconvolution input.
#'   For single-sample: named list of method matrices.
#'   For multi-sample: nested list by sample.
#' @param markers Optional marker list. If provided with an existing AEGIS object,
#'   it replaces `x$markers`.
#' @param do_marker Logical; run `audit_marker()`.
#' @param do_spatial Logical; run `audit_spatial()`.
#' @param do_compare Logical; run `compare_methods()`.
#' @param do_consensus Logical; run `compute_consensus()`.
#' @param strict Logical; strict validation for multi-sample constructor path.
#'
#' @return An `aegis` or `aegis_multi` object.
#' @export
run_aegis <- function(
    x,
    deconv = NULL,
    markers = NULL,
    do_marker = TRUE,
    do_spatial = TRUE,
    do_compare = TRUE,
    do_consensus = TRUE,
    strict = TRUE) {
  flags <- list(
    do_marker = do_marker,
    do_spatial = do_spatial,
    do_compare = do_compare,
    do_consensus = do_consensus,
    strict = strict
  )
  bad_flags <- names(flags)[!vapply(flags, function(v) is.logical(v) && length(v) == 1L && !is.na(v), logical(1))]
  if (length(bad_flags) > 0L) {
    stop(
      sprintf("These arguments must be single TRUE/FALSE values: %s", paste(bad_flags, collapse = ", ")),
      call. = FALSE
    )
  }

  obj <- NULL
  if (inherits(x, "aegis") || inherits(x, "aegis_multi")) {
    obj <- x
    if (!is.null(deconv)) {
      warning("`deconv` is ignored when `x` is already an AEGIS object.", call. = FALSE)
    }
    if (!is.null(markers)) {
      obj$markers <- markers
    }
  } else if (inherits(x, "Seurat")) {
    if (is.null(deconv)) {
      stop("`deconv` is required when `x` is a Seurat object.", call. = FALSE)
    }
    obj <- as_aegis(x, deconv = deconv, markers = markers)
  } else if (is.list(x) && length(x) > 0L && all(vapply(x, inherits, logical(1), what = "Seurat"))) {
    if (is.null(deconv)) {
      stop("`deconv` is required when `x` is a list of Seurat objects.", call. = FALSE)
    }
    obj <- as_aegis_multi(x, deconv = deconv, markers = markers, strict = strict)
  } else {
    stop("`x` must be a Seurat object, list of Seurat objects, aegis, or aegis_multi.", call. = FALSE)
  }

  obj <- audit_basic(obj)
  if (isTRUE(do_marker)) {
    obj <- audit_marker(obj)
  }
  if (isTRUE(do_spatial)) {
    obj <- audit_spatial(obj)
  }
  if (isTRUE(do_compare)) {
    obj <- compare_methods(obj)
  }
  if (isTRUE(do_consensus)) {
    obj <- compute_consensus(obj)
  }

  obj
}
