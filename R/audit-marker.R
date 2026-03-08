#' Run marker-based concordance audit
#'
#' Computes marker expression scores per cell type and compares them with
#' deconvolution abundance estimates.
#'
#' @param x An `aegis` object.
#' @param markers Named list of marker genes by cell type. Defaults to `x$markers`.
#' @param assay Assay name to use from Seurat.
#'
#' @return The modified `aegis` object with `x$audit$marker` populated.
#' @export
audit_marker <- function(x, markers = x$markers, assay = "Spatial") {
  assert_is_aegis(x)

  if (is.null(markers) || !is.list(markers) || length(markers) == 0L) {
    warning("No marker list provided. Skipping marker audit.", call. = FALSE)
    x$audit$marker <- list(
      marker_scores = NULL,
      concordance = data.frame(),
      summary_table = data.frame(),
      marker_usage = data.frame()
    )
    return(x)
  }

  if (!(assay %in% names(x$seu@assays))) {
    stop(sprintf("Assay '%s' not found in Seurat object.", assay), call. = FALSE)
  }

  expr <- tryCatch(
    Seurat::GetAssayData(x$seu, assay = assay, layer = "data"),
    error = function(e) Seurat::GetAssayData(x$seu, assay = assay, slot = "data")
  )
  if (nrow(expr) == 0L) {
    counts <- tryCatch(
      Seurat::GetAssayData(x$seu, assay = assay, layer = "counts"),
      error = function(e) Seurat::GetAssayData(x$seu, assay = assay, slot = "counts")
    )
    expr <- log1p(t(t(counts) / pmax(1, Matrix::colSums(counts))) * 1e4)
  }

  marker_scores <- list()
  marker_usage <- lapply(names(markers), function(ct) {
    genes <- unique(markers[[ct]])
    available <- intersect(genes, rownames(expr))

    if (length(available) == 0L) {
      marker_scores[[ct]] <<- rep(NA_real_, ncol(expr))
    } else {
      marker_scores[[ct]] <<- Matrix::colMeans(expr[available, , drop = FALSE])
    }

    data.frame(
      celltype = ct,
      n_markers_requested = length(genes),
      n_markers_used = length(available),
      stringsAsFactors = FALSE
    )
  })

  marker_score_mat <- do.call(cbind, marker_scores)
  colnames(marker_score_mat) <- names(markers)
  rownames(marker_score_mat) <- colnames(x$seu)

  usage_tbl <- dplyr::bind_rows(marker_usage)

  concordance <- dplyr::bind_rows(lapply(names(x$deconv), function(method) {
    mat <- x$deconv[[method]]
    overlap_ct <- intersect(colnames(mat), colnames(marker_score_mat))

    dplyr::bind_rows(lapply(overlap_ct, function(ct) {
      score <- marker_score_mat[, ct]
      pred <- mat[, ct]
      ok <- is.finite(score) & is.finite(pred)

      cor_val <- if (sum(ok) >= 3L) {
        stats::cor(score[ok], pred[ok], method = "pearson")
      } else {
        NA_real_
      }

      used_n <- usage_tbl$n_markers_used[match(ct, usage_tbl$celltype)]
      data.frame(
        method = method,
        celltype = ct,
        correlation = cor_val,
        n_markers_used = used_n,
        stringsAsFactors = FALSE
      )
    }))
  }))

  summary_table <- concordance |>
    dplyr::group_by(.data$method) |>
    dplyr::summarise(
      mean_correlation = mean(.data$correlation, na.rm = TRUE),
      median_correlation = stats::median(.data$correlation, na.rm = TRUE),
      n_celltypes = dplyr::n(),
      .groups = "drop"
    )

  x$audit$marker <- list(
    marker_scores = marker_score_mat,
    concordance = concordance,
    summary_table = summary_table,
    marker_usage = usage_tbl,
    assay = assay
  )

  if (any(usage_tbl$n_markers_used == 0L)) {
    missing_ct <- usage_tbl$celltype[usage_tbl$n_markers_used == 0L]
    warning(
      sprintf("Some marker sets had no genes present and were skipped: %s", paste(missing_ct, collapse = ", ")),
      call. = FALSE
    )
  }

  x
}
