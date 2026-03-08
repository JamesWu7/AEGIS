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
  assert_is_seurat(x$seu, "x$seu")

  if (is.null(x$deconv) || !is.list(x$deconv) || length(x$deconv) == 0L) {
    stop("`x$deconv` must be a non-empty list.", call. = FALSE)
  }

  if (is.null(markers)) {
    warning("No marker list provided. Skipping marker audit.", call. = FALSE)
    x$audit$marker <- list(
      marker_scores = NULL,
      concordance = data.frame(),
      marker_usage = data.frame(),
      summary = data.frame(),
      summary_table = data.frame(),
      assay = assay,
      status = "no_markers"
    )
    return(x)
  }
  if (!is.list(markers)) {
    stop("`markers` must be NULL or a named list.", call. = FALSE)
  }
  marker_names <- names(markers)
  if (is.null(marker_names) || any(trimws(marker_names) == "")) {
    stop("`markers` must be a named list with non-empty cell-type names.", call. = FALSE)
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
    lib_size <- Matrix::colSums(counts)
    lib_size[lib_size <= 0] <- 1
    scale_diag <- Matrix::Diagonal(x = 1e4 / lib_size)
    expr <- log1p(counts %*% scale_diag)
    rownames(expr) <- rownames(counts)
    colnames(expr) <- colnames(counts)
  }
  if (nrow(expr) == 0L || ncol(expr) == 0L) {
    stop(sprintf("Assay '%s' has no usable expression data.", assay), call. = FALSE)
  }

  spots <- colnames(x$seu)
  expr <- expr[, spots, drop = FALSE]

  marker_scores <- list()
  marker_usage <- lapply(names(markers), function(ct) {
    genes <- unique(as.character(markers[[ct]]))
    genes <- genes[!is.na(genes) & trimws(genes) != ""]
    available <- intersect(genes, rownames(expr))

    if (length(available) == 0L) {
      marker_scores[[ct]] <<- rep(NA_real_, length(spots))
    } else {
      marker_scores[[ct]] <<- Matrix::colMeans(expr[available, , drop = FALSE])
    }

    data.frame(
      celltype = ct,
      n_markers_requested = length(genes),
      n_markers_used = length(available),
      markers_used = paste(available, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  marker_score_mat <- do.call(cbind, marker_scores)
  colnames(marker_score_mat) <- names(markers)
  rownames(marker_score_mat) <- spots

  usage_tbl <- dplyr::bind_rows(marker_usage)

  concordance <- dplyr::bind_rows(lapply(names(x$deconv), function(method) {
    mat <- x$deconv[[method]]
    if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
      stop(sprintf("Method '%s' must be a non-empty numeric matrix.", method), call. = FALSE)
    }
    if (!identical(rownames(mat), spots)) {
      stop(sprintf("Method '%s' rownames are not aligned to Seurat spots.", method), call. = FALSE)
    }
    overlap_ct <- intersect(colnames(mat), colnames(marker_score_mat))

    if (length(overlap_ct) == 0L) {
      return(data.frame())
    }

    dplyr::bind_rows(lapply(overlap_ct, function(ct) {
      score <- marker_score_mat[, ct]
      pred <- mat[, ct]
      ok <- is.finite(score) & is.finite(pred)

      pearson_cor <- if (sum(ok) >= 3L) {
        stats::cor(score[ok], pred[ok], method = "pearson")
      } else {
        NA_real_
      }
      spearman_cor <- if (sum(ok) >= 3L) {
        stats::cor(score[ok], pred[ok], method = "spearman")
      } else {
        NA_real_
      }

      used_n <- usage_tbl$n_markers_used[match(ct, usage_tbl$celltype)]
      data.frame(
        method = method,
        celltype = ct,
        pearson_cor = pearson_cor,
        spearman_cor = spearman_cor,
        n_markers_used = used_n,
        n_spots_used = sum(ok),
        stringsAsFactors = FALSE
      )
    }))
  }))

  if (nrow(concordance) > 0L) {
    summary_tbl <- concordance |>
      dplyr::group_by(.data$method) |>
      dplyr::summarise(
        mean_pearson = mean(.data$pearson_cor, na.rm = TRUE),
        mean_spearman = mean(.data$spearman_cor, na.rm = TRUE),
        n_celltypes = dplyr::n(),
        .groups = "drop"
      )
  } else {
    summary_tbl <- data.frame()
  }

  x$audit$marker <- list(
    marker_scores = marker_score_mat,
    concordance = concordance,
    marker_usage = usage_tbl,
    summary = summary_tbl,
    summary_table = summary_tbl,
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
