#' Plot audit outputs
#'
#' @param x An `aegis` object.
#' @param type One of `sumdev`, `dominance`, `entropy`, `marker`, `smoothness`.
#'
#' @return A ggplot object.
#' @export
plot_audit <- function(x, type = c("sumdev", "dominance", "entropy", "marker", "smoothness")) {
  assert_is_aegis(x)
  type <- match.arg(type)

  if (type %in% c("sumdev", "dominance", "entropy")) {
    if (is.null(x$audit$basic$spot_metrics)) {
      stop("Basic audit not found. Run audit_basic() first.", call. = FALSE)
    }
    coords <- get_spatial_coordinates(x$seu)
    dat <- dplyr::left_join(x$audit$basic$spot_metrics, coords, by = c("spot" = "spot"))

    metric <- switch(type,
      sumdev = "sum_dev",
      dominance = "dominance",
      entropy = "entropy"
    )

    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$imagecol, y = .data$imagerow, color = .data[[metric]])) +
      ggplot2::geom_point(size = 0.8) +
      ggplot2::facet_wrap(~method) +
      ggplot2::scale_y_reverse() +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(x = "imagecol", y = "imagerow", color = metric, title = sprintf("%s spatial map", metric))

    return(p)
  }

  if (type == "marker") {
    if (is.null(x$audit$marker$concordance)) {
      stop("Marker audit not found. Run audit_marker() first.", call. = FALSE)
    }

    p <- ggplot2::ggplot(x$audit$marker$concordance,
      ggplot2::aes(x = .data$celltype, y = .data$correlation, fill = .data$method)
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "Marker concordance by method", x = "Cell type", y = "Pearson correlation")

    return(p)
  }

  if (type == "smoothness") {
    if (is.null(x$audit$spatial$method_summary)) {
      stop("Spatial audit not found. Run audit_spatial() first.", call. = FALSE)
    }

    p <- ggplot2::ggplot(x$audit$spatial$method_summary,
      ggplot2::aes(x = .data$celltype, y = .data$mean_local_inconsistency, fill = .data$method)
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "Spatial local inconsistency", x = "Cell type", y = "Mean local inconsistency")

    return(p)
  }

  stop("Unsupported `type`.", call. = FALSE)
}
