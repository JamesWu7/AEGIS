#' Compare deconvolution methods
#'
#' Computes pairwise method agreement by cell type and spot-level agreement
#' across methods.
#'
#' @param x An `aegis` object.
#'
#' @return Modified `aegis` object with `x$consensus$comparison` populated.
#' @export
compare_methods <- function(x) {
  assert_is_aegis(x)
  if (is.null(x$deconv) || !is.list(x$deconv) || length(x$deconv) == 0L) {
    stop("`x$deconv` must be a non-empty list.", call. = FALSE)
  }

  methods <- names(x$deconv)
  if (length(methods) < 2L) {
    stop("At least two methods are required for comparison.", call. = FALSE)
  }
  spots <- colnames(x$seu)

  for (method in methods) {
    mat <- x$deconv[[method]]
    if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
      stop(sprintf("Method '%s' must be a non-empty numeric matrix.", method), call. = FALSE)
    }
    if (!identical(rownames(mat), spots)) {
      stop(sprintf("Method '%s' rownames are not aligned to Seurat spots.", method), call. = FALSE)
    }
  }

  pairs <- utils::combn(methods, 2, simplify = FALSE)

  pairwise_celltype_cor <- dplyr::bind_rows(lapply(pairs, function(p) {
    m1 <- p[[1]]
    m2 <- p[[2]]
    mat1 <- x$deconv[[m1]]
    mat2 <- x$deconv[[m2]]

    shared_ct <- intersect(colnames(mat1), colnames(mat2))
    if (length(shared_ct) == 0L) {
      return(data.frame())
    }

    dplyr::bind_rows(lapply(shared_ct, function(ct) {
      cor_val <- stats::cor(mat1[, ct], mat2[, ct], use = "pairwise.complete.obs")
      data.frame(
        method_1 = m1,
        method_2 = m2,
        celltype = ct,
        correlation = cor_val,
        n_spots = nrow(mat1),
        stringsAsFactors = FALSE
      )
    }))
  }))

  if (nrow(pairwise_celltype_cor) == 0L) {
    stop("No overlapping cell types found across methods.", call. = FALSE)
  }
  pairwise_celltype_cor <- tidyr::drop_na(pairwise_celltype_cor, "correlation")

  pairwise_summary <- pairwise_celltype_cor |>
    dplyr::group_by(.data$method_1, .data$method_2) |>
    dplyr::summarise(
      mean_correlation = mean(.data$correlation, na.rm = TRUE),
      n_shared_celltypes = dplyr::n(),
      .groups = "drop"
    )

  heatmap_matrix <- matrix(1, nrow = length(methods), ncol = length(methods),
    dimnames = list(methods, methods)
  )
  for (i in seq_len(nrow(pairwise_summary))) {
    m1 <- pairwise_summary$method_1[[i]]
    m2 <- pairwise_summary$method_2[[i]]
    v <- pairwise_summary$mean_correlation[[i]]
    heatmap_matrix[m1, m2] <- v
    heatmap_matrix[m2, m1] <- v
  }

  celltype_agreement <- pairwise_celltype_cor |>
    dplyr::group_by(.data$celltype) |>
    dplyr::summarise(
      average_agreement = mean(.data$correlation, na.rm = TRUE),
      median_agreement = stats::median(.data$correlation, na.rm = TRUE),
      n_pairs = dplyr::n(),
      .groups = "drop"
    )

  # Spot-level agreement from cross-method mean absolute deviation over cell types.
  all_ct <- sort(unique(unlist(lapply(x$deconv, colnames))))

  spot_agreement <- data.frame(
    spot = spots,
    mean_mad = NA_real_,
    agreement = NA_real_,
    stringsAsFactors = FALSE
  )

  mean_mad <- vapply(spots, function(s) {
    ct_dev <- vapply(all_ct, function(ct) {
      vals <- vapply(x$deconv, function(mat) {
        if (ct %in% colnames(mat)) mat[s, ct] else NA_real_
      }, numeric(1))
      vals <- vals[is.finite(vals)]
      if (length(vals) < 2L) return(NA_real_)
      mean(abs(vals - mean(vals)))
    }, numeric(1))
    mean(ct_dev, na.rm = TRUE)
  }, numeric(1))

  spot_agreement$mean_mad <- mean_mad
  spot_agreement$agreement <- 1 / (1 + mean_mad)

  compared_celltypes <- dplyr::bind_rows(lapply(pairs, function(p) {
    m1 <- p[[1]]
    m2 <- p[[2]]
    shared_ct <- intersect(colnames(x$deconv[[m1]]), colnames(x$deconv[[m2]]))
    data.frame(
      method_1 = m1,
      method_2 = m2,
      shared_celltypes = paste(shared_ct, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))

  x$consensus$comparison <- list(
    pairwise_celltype_cor = pairwise_celltype_cor,
    pairwise_celltype = pairwise_celltype_cor,
    pairwise_summary = pairwise_summary,
    heatmap_matrix = heatmap_matrix,
    celltype_agreement = celltype_agreement,
    spot_agreement = spot_agreement,
    compared_celltypes = compared_celltypes,
    summary = pairwise_summary
  )

  x
}
