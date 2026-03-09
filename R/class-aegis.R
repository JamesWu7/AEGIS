#' AEGIS analysis object
#'
#' `aegis` is a lightweight S3 container around a Seurat object and analysis state.
#'
#' @param x An object.
#'
#' @return Logical scalar.
#' @export
is_aegis <- function(x) {
  inherits(x, "aegis")
}

#' Check if object is an AEGIS multi-sample object
#'
#' @param x An object.
#'
#' @return Logical scalar.
#' @export
is_aegis_multi <- function(x) {
  inherits(x, "aegis_multi")
}

#' @export
print.aegis <- function(x, ...) {
  methods_n <- length(x$deconv)
  spots_n <- ncol(x$seu)
  genes_n <- nrow(x$seu)
  audits <- names(x$audit)
  audits_label <- if (length(audits) == 0L) "none" else paste(audits, collapse = ", ")

  cat("<aegis>\n")
  cat(sprintf("  Spots: %d\n", spots_n))
  cat(sprintf("  Genes: %d\n", genes_n))
  cat(sprintf("  Methods: %d (%s)\n", methods_n, paste(names(x$deconv), collapse = ", ")))
  cat(sprintf("  Marker sets: %s\n", if (is.null(x$markers)) "none" else length(x$markers)))
  cat(sprintf("  Audits: %s\n", audits_label))
  invisible(x)
}

#' @export
print.aegis_multi <- function(x, ...) {
  sample_ids <- names(x$seu_list %||% list())
  n_samples <- length(sample_ids)
  methods <- unique(unlist(lapply(x$deconv %||% list(), names)))
  methods_label <- if (length(methods) == 0L) "none" else paste(methods, collapse = ", ")
  audits <- names(x$audit)
  audits_label <- if (length(audits) == 0L) "none" else paste(audits, collapse = ", ")

  cat("<aegis_multi>\n")
  cat(sprintf("  Samples: %d (%s)\n", n_samples, paste(sample_ids, collapse = ", ")))
  cat(sprintf("  Methods: %s\n", methods_label))
  cat(sprintf("  Marker sets: %s\n", if (is.null(x$markers)) "none" else length(x$markers)))
  cat(sprintf("  Audits: %s\n", audits_label))
  invisible(x)
}
