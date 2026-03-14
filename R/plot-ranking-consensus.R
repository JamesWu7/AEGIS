#' Plot method ranking summary
#'
#' @param x An `aegis` object.
#' @param sample Optional sample ID for multi-sample objects.
#' @param palette Palette family: `nature`, `viridis`, `scico`, or `brewer`.
#' @param base_size Base font size.
#'
#' @return A ggplot object.
#' @export
plot_method_ranking <- function(
    x,
    sample = NULL,
    palette = "nature",
    base_size = 12) {
  x <- resolve_plot_input_aegis(x, sample = sample, fun_name = "plot_method_ranking")
  assert_is_aegis(x)
  tbl <- get_method_ranking_table(x)
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) {
    # Keep plotting ergonomic: attempt to derive ranking from available evidence.
    x <- tryCatch(
      {
        x2 <- x
        if (is.null(get_method_evidence_table(x2))) {
          x2 <- score_methods(x2)
        }
        rank_methods(x2, method = "mean_rank")
      },
      error = function(e) {
        stop(
          sprintf("Method ranking not found and auto-ranking failed: %s", conditionMessage(e)),
          call. = FALSE
        )
      }
    )
    tbl <- get_method_ranking_table(x)
    if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) {
      stop("Method ranking not found. Run rank_methods() first.", call. = FALSE)
    }
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
  dat$rank_label <- paste0("#", dat$overall_rank)
  n_methods <- nrow(dat)
  score_min <- min(dat$overall_score, na.rm = TRUE)
  score_min <- if (is.finite(score_min)) min(score_min, 0) else 0
  score_span <- diff(range(dat$overall_score, na.rm = TRUE))
  if (!is.finite(score_span) || score_span == 0) score_span <- 1
  show_rank_label <- n_methods <= 16L
  point_size <- if (n_methods > 12L) 3.0 else 3.6
  line_width <- if (n_methods > 12L) 0.86 else 1.0
  label_size <- if (n_methods > 12L) base_size * 0.20 else base_size * 0.23

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$overall_score, y = .data$method, color = .data$recommendation)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = score_min, xend = .data$overall_score, y = .data$method, yend = .data$method),
      linewidth = line_width,
      alpha = 0.75
    ) +
    ggplot2::geom_point(size = point_size, alpha = 0.96) +
    scale_color_aegis(palette = palette, type = "categorical") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.14))) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_aegis(base_size = base_size) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::labs(
      title = "Method Ranking",
      subtitle = "Aggregated from marker, spatial, agreement, and stability evidence",
      x = "Aggregated score",
      y = "Method",
      color = "Recommendation"
    )

  if (isTRUE(show_rank_label)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$rank_label),
      nudge_x = 0.02 * score_span,
      hjust = 0,
      size = label_size,
      color = "#222222",
      show.legend = FALSE
    )
  }

  p
}

#' Plot spot-level consensus disagreement map
#'
#' @param x An `aegis` object.
#' @param sample Optional sample ID for multi-sample objects.
#' @param palette Palette family: `nature`, `viridis`, `scico`, or `brewer`.
#' @param base_size Base font size.
#'
#' @return A plot object (Seurat spatial plot when available, otherwise ggplot).
#' @export
plot_disagreement_map <- function(
    x,
    sample = NULL,
    palette = "nature",
    base_size = 12) {
  x <- resolve_plot_input_aegis(x, sample = sample, fun_name = "plot_disagreement_map")
  assert_is_aegis(x)
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
#' @param sample Optional sample ID for multi-sample objects.
#' @param palette Palette family: `nature`, `viridis`, `scico`, or `brewer`.
#' @param base_size Base font size.
#'
#' @return A plot object (Seurat spatial plot when available, otherwise ggplot).
#' @export
plot_consensus_confidence <- function(
    x,
    sample = NULL,
    palette = "nature",
    base_size = 12) {
  x <- resolve_plot_input_aegis(x, sample = sample, fun_name = "plot_consensus_confidence")
  assert_is_aegis(x)
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
    pal_vals <- get_plot_palette(palette = palette, type = "continuous", n = 256L)
    col_low <- pal_vals[[1L]]
    col_high <- pal_vals[[length(pal_vals)]]
    p <- tryCatch(
      Seurat::SpatialFeaturePlot(
        seu,
        features = metric_name,
        cols = c(col_low, col_high),
        image.alpha = 0.25,
        pt.size.factor = 1.85
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
  point_size <- auto_spatial_point_size(nrow(spatial_dat))
  ggplot2::ggplot(spatial_dat, ggplot2::aes(x = .data$x, y = .data$y, color = .data$value)) +
    ggplot2::geom_point(size = point_size, alpha = 0.95) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    scale_color_aegis(palette = palette, type = "continuous") +
    theme_aegis_spatial(base_size = base_size) +
    ggplot2::labs(title = title, subtitle = subtitle, color = legend_title)
}
