#' Compute consensus deconvolution profile
#'
#' Aggregates deconvolution outputs across methods (mean aggregation for MVP)
#' and computes disagreement/stability summaries.
#'
#' @param x An `aegis` object.
#'
#' @return Modified `aegis` object with `x$consensus$result` populated.
#' @export
compute_consensus <- function(x) {
  if (is_multi_sample_context(x)) {
    sample_results <- iterate_aegis_samples(x, compute_consensus)
    by_sample <- lapply(sample_results, function(obj) obj$consensus$result)

    spot_confidence <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]$spot_confidence %||% data.frame()
      if (nrow(tbl) == 0L) return(tbl)
      tbl$sample_id <- sid
      tbl
    }))
    celltype_stability <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]$celltype_stability %||% data.frame()
      if (nrow(tbl) == 0L) return(tbl)
      tbl$sample_id <- sid
      tbl
    }))
    shared_celltypes <- lapply(by_sample, function(z) z$shared_celltypes %||% character())

    x$consensus$result <- list(
      by_sample = by_sample,
      spot_confidence = spot_confidence,
      celltype_stability = celltype_stability,
      shared_celltypes = shared_celltypes
    )
    return(x)
  }

  assert_is_aegis(x)
  if (is.null(x$deconv) || !is.list(x$deconv) || length(x$deconv) == 0L) {
    stop("`x$deconv` must be a non-empty list.", call. = FALSE)
  }

  methods <- names(x$deconv)
  spots <- colnames(x$seu)
  if (is.null(spots) || length(spots) == 0L) {
    stop("Seurat object has no spot names.", call. = FALSE)
  }
  for (method in methods) {
    mat <- x$deconv[[method]]
    if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
      stop(sprintf("Method '%s' must be a non-empty numeric matrix.", method), call. = FALSE)
    }
    if (!identical(rownames(mat), spots)) {
      stop(sprintf("Method '%s' rownames are not aligned to Seurat spots.", method), call. = FALSE)
    }
  }

  shared_celltypes <- Reduce(intersect, lapply(x$deconv, colnames))
  if (length(shared_celltypes) == 0L) {
    stop("No shared cell types found across methods for consensus.", call. = FALSE)
  }
  shared_celltypes <- sort(shared_celltypes)

  values_by_method <- lapply(methods, function(method) {
    x$deconv[[method]][spots, shared_celltypes, drop = FALSE]
  })
  names(values_by_method) <- methods
  arr <- simplify2array(values_by_method)

  consensus <- apply(arr, c(1, 2), mean)
  disagreement <- apply(arr, c(1, 2), stats::sd)

  rownames(consensus) <- spots
  colnames(consensus) <- shared_celltypes
  rownames(disagreement) <- spots
  colnames(disagreement) <- shared_celltypes

  consensus[!is.finite(consensus)] <- NA_real_
  disagreement[!is.finite(disagreement)] <- NA_real_

  rs <- rowSums(consensus)
  rs[!is.finite(rs) | rs <= 0] <- 1
  consensus <- consensus / rs

  mean_sd <- rowMeans(disagreement, na.rm = TRUE)
  mean_sd[!is.finite(mean_sd)] <- NA_real_
  spot_confidence <- data.frame(
    spot = spots,
    mean_disagreement = mean_sd,
    confidence = 1 / (1 + mean_sd),
    stringsAsFactors = FALSE
  )

  ct_sd <- colMeans(disagreement, na.rm = TRUE)
  ct_sd[!is.finite(ct_sd)] <- NA_real_
  celltype_stability <- data.frame(
    celltype = shared_celltypes,
    mean_disagreement = ct_sd,
    stability = 1 / (1 + ct_sd),
    stringsAsFactors = FALSE
  )

  x$consensus$result <- list(
    consensus_matrix = consensus,
    method_disagreement = disagreement,
    spot_confidence = spot_confidence,
    celltype_stability = celltype_stability,
    shared_celltypes = shared_celltypes,
    methods = methods
  )

  x
}
