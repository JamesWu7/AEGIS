#' Run basic quality audits on deconvolution matrices
#'
#' Computes per-spot diagnostics and method-level sparsity summaries.
#'
#' @param x An `aegis` object.
#' @param threshold Numeric threshold for counting detected cell types.
#'
#' @return The modified `aegis` object with `x$audit$basic` populated.
#' @export
audit_basic <- function(x, threshold = 0.05) {
  assert_is_aegis(x)
  if (is.null(x$deconv) || !is.list(x$deconv) || length(x$deconv) == 0L) {
    stop("`x$deconv` must be a non-empty list of deconvolution matrices.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) || threshold < 0) {
    stop("`threshold` must be a single non-negative numeric value.", call. = FALSE)
  }

  per_method <- lapply(names(x$deconv), function(method) {
    mat <- x$deconv[[method]]
    if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
      stop(sprintf("Method '%s' must be a non-empty numeric matrix.", method), call. = FALSE)
    }
    if (is.null(rownames(mat)) || any(trimws(rownames(mat)) == "")) {
      stop(sprintf("Method '%s' must contain non-empty row names.", method), call. = FALSE)
    }
    if (anyNA(mat)) {
      stop(sprintf("Method '%s' contains NA values.", method), call. = FALSE)
    }

    sum_dev <- abs(rowSums(mat) - 1)
    dominance <- apply(mat, 1, max)
    entropy <- row_entropy(mat)
    n_detected <- rowSums(mat > threshold)

    data.frame(
      spot = rownames(mat),
      method = method,
      sum_dev = sum_dev,
      dominance = dominance,
      entropy = entropy,
      n_detected_types = n_detected,
      stringsAsFactors = FALSE
    )
  })
  names(per_method) <- names(x$deconv)

  summary <- dplyr::bind_rows(lapply(names(x$deconv), function(method) {
    mat <- x$deconv[[method]]
    spot_tbl <- per_method[[method]]
    data.frame(
      method = method,
      n_spots = nrow(mat),
      n_celltypes = ncol(mat),
      zero_fraction = mean(mat == 0),
      near_zero_fraction = mean(mat <= threshold),
      mean_dominance = mean(spot_tbl$dominance),
      mean_entropy = mean(spot_tbl$entropy),
      mean_n_detected_types = mean(spot_tbl$n_detected_types),
      mean_sum_dev = mean(spot_tbl$sum_dev),
      stringsAsFactors = FALSE
    )
  }))

  x$audit$basic <- list(
    per_method_spot_metrics = per_method,
    spot_metrics = dplyr::bind_rows(per_method),
    summary = summary,
    summary_table = summary,
    sparsity_summary = summary[, c("method", "zero_fraction", "near_zero_fraction"), drop = FALSE],
    threshold = threshold
  )

  x
}
