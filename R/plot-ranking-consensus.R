#' Plot method ranking summary
#'
#' @param x An `aegis` object.
#' @param palette Palette family: `nature`, `viridis`, `scico`, or `brewer`.
#' @param base_size Base font size.
#'
#' @return A ggplot object.
#' @export
plot_method_ranking <- function(
    x,
    palette = "nature",
    base_size = 12) {
  assert_is_aegis(x)
  tbl <- get_method_ranking_table(x)
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) {
    stop("Method ranking not found. Run rank_methods() first.", call. = FALSE)
  }
  if (!all(c("method", "overall_rank") %in% colnames(tbl))) {
    stop("Method ranking table must include `method` and `overall_rank`.", call. = FALSE)
  }

  dat <- tbl
  if (!("overall_score" %in% colnames(dat))) {
    dat$overall_score <- 1 / dat$overall_rank
  }
  if (!("recommendation" %in% colnames(dat))) {
    dat$recommendation <- recommendation_from_rank(dat$overall_rank)
  }

  dat <- dat[order(dat$overall_rank), , drop = FALSE]
  dat$method <- factor(dat$method, levels = rev(dat$method))

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$overall_score, y = .data$method, color = .data$recommendation)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = min(.data$overall_score, na.rm = TRUE), xend = .data$overall_score, y = .data$method, yend = .data$method),
      linewidth = 0.8,
      alpha = 0.7
    ) +
    ggplot2::geom_point(size = 3.2, alpha = 0.95) +
    scale_color_aegis(palette = palette, type = "categorical") +
    theme_aegis(base_size = base_size) +
    ggplot2::labs(
      title = "Method Ranking",
      subtitle = "Aggregated from marker, spatial, agreement, and stability evidence",
      x = "Aggregated score",
      y = "Method",
      color = "Recommendation"
    )
}

#' Plot spot-level consensus disagreement map
#'
#' @param x An `aegis` object.
#' @param palette Palette family: `nature`, `viridis`, `scico`, or `brewer`.
#' @param base_size Base font size.
#'
#' @return A plot object (Seurat spatial plot when available, otherwise ggplot).
#' @export
plot_disagreement_map <- function(
    x,
    palette = "nature",
    base_size = 12) {
  assert_is_aegis(x)
  if (is_multi_sample_context(x)) {
    stop("Multi-sample object detected. Use split_aegis_by_sample() and plot one sample at a time.", call. = FALSE)
  }
  if (is.null(x$consensus$result$method_disagreement)) {
    stop("Consensus disagreement not found. Run compute_consensus() first.", call. = FALSE)
  }

  dis <- x$consensus$result$method_disagreement
  if (!is.matrix(dis) || nrow(dis) == 0L || ncol(dis) == 0L || is.null(rownames(dis))) {
    stop("`x$consensus$result$method_disagreement` must be a non-empty matrix with spot rownames.", call. = FALSE)
  }
  spot_disagreement <- rowMeans(dis, na.rm = TRUE)
  spot_disagreement[!is.finite(spot_disagreement)] <- NA_real_
  names(spot_disagreement) <- rownames(dis)

  plot_spatial_metric(
    seu = x$seu,
    spot_values = spot_disagreement,
    metric_name = ".aegis_method_disagreement",
    legend_title = "Disagreement",
    title = "Method Disagreement Map",
    subtitle = "Higher values indicate stronger cross-method disagreement",
    palette = palette,
    base_size = base_size
  )
}

#' Plot spot-level consensus confidence map
#'
#' @param x An `aegis` object.
#' @param palette Palette family: `nature`, `viridis`, `scico`, or `brewer`.
#' @param base_size Base font size.
#'
#' @return A plot object (Seurat spatial plot when available, otherwise ggplot).
#' @export
plot_consensus_confidence <- function(
    x,
    palette = "nature",
    base_size = 12) {
  assert_is_aegis(x)
  if (is_multi_sample_context(x)) {
    stop("Multi-sample object detected. Use split_aegis_by_sample() and plot one sample at a time.", call. = FALSE)
  }
  if (is.null(x$consensus$result$spot_confidence)) {
    stop("Consensus confidence not found. Run compute_consensus() first.", call. = FALSE)
  }

  conf <- x$consensus$result$spot_confidence
  if (!is.data.frame(conf) || !all(c("spot", "confidence") %in% colnames(conf))) {
    stop("`x$consensus$result$spot_confidence` must include `spot` and `confidence`.", call. = FALSE)
  }

  vals <- conf$confidence
  vals[!is.finite(vals)] <- NA_real_
  names(vals) <- conf$spot

  plot_spatial_metric(
    seu = x$seu,
    spot_values = vals,
    metric_name = ".aegis_consensus_confidence",
    legend_title = "Confidence",
    title = "Consensus Confidence Map",
    subtitle = "Higher values indicate more reliable integrated predictions",
    palette = palette,
    base_size = base_size
  )
}

#' @keywords internal
plot_spatial_metric <- function(
    seu,
    spot_values,
    metric_name,
    legend_title,
    title,
    subtitle,
    palette = "nature",
    base_size = 12) {
  assert_is_seurat(seu)
  spots <- colnames(seu)
  vals <- as.numeric(spot_values[match(spots, names(spot_values))])

  seu[[metric_name]] <- vals
  has_image <- FALSE
  has_image <- tryCatch(length(Seurat::Images(seu)) > 0L, error = function(e) FALSE)

  if (isTRUE(has_image)) {
    p <- tryCatch(
      Seurat::SpatialFeaturePlot(
        seu,
        features = metric_name,
        cols = get_plot_palette(palette = palette, type = "continuous", n = 2L),
        image.alpha = 0.25,
        pt.size.factor = 1.6
      ),
      error = function(e) NULL
    )
    if (!is.null(p)) {
      p <- p + ggplot2::labs(
        title = title,
        subtitle = subtitle,
        fill = legend_title
      ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = base_size * 1.2, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = base_size * 0.95)
        )
      return(p)
    }
  }

  dat <- data.frame(
    spot = spots,
    value = vals,
    stringsAsFactors = FALSE
  )
  spatial_dat <- assemble_spatial_plot_data(seu, dat)
  ggplot2::ggplot(spatial_dat, ggplot2::aes(x = .data$x, y = .data$y, color = .data$value)) +
    ggplot2::geom_point(size = if (nrow(spatial_dat) > 5000L) 0.45 else 0.7, alpha = 0.95) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    scale_color_aegis(palette = palette, type = "continuous") +
    theme_aegis_spatial(base_size = base_size) +
    ggplot2::labs(title = title, subtitle = subtitle, color = legend_title)
}
