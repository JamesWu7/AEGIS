#' Run basic quality audits on deconvolution matrices
#'
#' Computes per-spot diagnostics and method-level sparsity summaries.
#'
#' @param x An `aegis` object.
#' @param detection_threshold Numeric threshold for counting detected cell types.
#'
#' @return The modified `aegis` object with `x$audit$basic` populated.
#' @export
audit_basic <- function(x, detection_threshold = 0.05) {
  assert_is_aegis(x)

  per_method <- lapply(names(x$deconv), function(method) {
    mat <- x$deconv[[method]]
    sum_dev <- abs(rowSums(mat) - 1)
    dominance <- apply(mat, 1, max)
    entropy <- row_entropy(mat)
    n_detected <- rowSums(mat > detection_threshold)

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

  spot_metrics <- dplyr::bind_rows(per_method)

  summary_table <- spot_metrics |>
    dplyr::group_by(.data$method) |>
    dplyr::summarise(
      mean_sum_dev = mean(.data$sum_dev, na.rm = TRUE),
      mean_dominance = mean(.data$dominance, na.rm = TRUE),
      mean_entropy = mean(.data$entropy, na.rm = TRUE),
      mean_n_detected = mean(.data$n_detected_types, na.rm = TRUE),
      .groups = "drop"
    )

  sparsity_summary <- dplyr::bind_rows(lapply(names(x$deconv), function(method) {
    mat <- x$deconv[[method]]
    data.frame(
      method = method,
      n_spots = nrow(mat),
      n_celltypes = ncol(mat),
      zero_fraction = mean(mat <= 1e-8),
      near_zero_fraction = mean(mat <= detection_threshold),
      stringsAsFactors = FALSE
    )
  }))

  x$audit$basic <- list(
    spot_metrics = spot_metrics,
    summary_table = summary_table,
    sparsity_summary = sparsity_summary,
    params = list(detection_threshold = detection_threshold)
  )

  x
}
