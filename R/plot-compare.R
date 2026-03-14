#' Plot method comparison and consensus outputs
#'
#' @param x An `aegis` object.
#' @param type One of `heatmap`, `spot_agreement`, `consensus_map`,
#'   `disagreement_map`, `confidence_map`, or `ranking`.
#' @param sample Optional sample ID for multi-sample objects.
#' @param palette Palette family: `nature`, `viridis`, `scico`, or `brewer`.
#' @param base_size Base font size for the plot theme.
#'
#' @return A ggplot object.
#' @export
plot_compare <- function(
    x,
    type = c("heatmap", "spot_agreement", "consensus_map", "disagreement_map", "confidence_map", "ranking"),
    sample = NULL,
    palette = "nature",
    base_size = 12) {
  x <- resolve_plot_input_aegis(x, sample = sample, fun_name = "plot_compare")

  assert_is_aegis(x)
  type <- match.arg(type)

  if (type == "heatmap") {
    if (is.null(x$consensus$comparison)) {
      stop("Comparison results not found. Run compare_methods() first.", call. = FALSE)
    }
    dat <- build_agreement_heatmap_data(x$consensus$comparison)
    dat <- dat[is.finite(dat$correlation), , drop = FALSE]
    if (nrow(dat) == 0L) {
      stop("No finite pairwise correlations available for heatmap plotting.", call. = FALSE)
    }

    pair_order <- dat |>
      dplyr::group_by(.data$method_pair) |>
      dplyr::summarise(mean_cor = mean(.data$correlation, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(.data$mean_cor)) |>
      dplyr::pull(.data$method_pair)
    n_pairs <- length(pair_order)
    n_methods <- length(unique(c(dat$method_1, dat$method_2)))

    pair_labels <- pair_order
    if (n_pairs > 12L) {
      pair_labels <- vapply(
        pair_order,
        function(lbl) gsub(" vs ", "\nvs\n", lbl, fixed = TRUE),
        character(1)
      )
    }

    dat$method_pair <- factor(dat$method_pair, levels = pair_order, labels = pair_labels)
    dat$celltype <- factor(dat$celltype, levels = sort(unique(dat$celltype)))

    axis_x_angle <- if (n_pairs > 24L) 90 else if (n_pairs > 12L) 60 else 30
    axis_x_size <- if (n_pairs > 36L) base_size * 0.56 else if (n_pairs > 24L) base_size * 0.62 else if (n_pairs > 12L) base_size * 0.72 else base_size * 0.82
    axis_y_size <- if (length(levels(dat$celltype)) > 10L) base_size * 0.72 else base_size * 0.84

    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$method_pair, y = .data$celltype, fill = .data$correlation)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_aegis(palette = palette, type = "diverging", limits = c(-1, 1)) +
      theme_aegis(base_size = base_size) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = axis_x_size, angle = axis_x_angle, hjust = 1, vjust = if (axis_x_angle == 90) 0.5 else 1),
        axis.text.y = ggplot2::element_text(size = axis_y_size),
        legend.key.width = grid::unit(if (n_pairs > 24L) 1.5 else 1.2, "lines")
      ) +
      ggplot2::labs(
        title = "Method-Pair Agreement by Cell Type",
        subtitle = sprintf("Correlation across shared spots (%d methods, %d pairs)", n_methods, n_pairs),
        x = "Method pair",
        y = "Cell type",
        fill = "Correlation"
      )
    if (length(unique(dat$celltype)) <= 8L && n_pairs <= 8L) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", .data$correlation)),
        size = base_size * 0.23,
        color = "#1F1F1F"
      )
    }
    return(p)
  }

  if (type == "spot_agreement") {
    if (is.null(x$consensus$comparison$spot_agreement)) {
      stop("Comparison spot agreement not found. Run compare_methods() first.", call. = FALSE)
    }
    dat <- x$consensus$comparison$spot_agreement
    if (!is.data.frame(dat) || !("spot" %in% colnames(dat))) {
      stop("`x$consensus$comparison$spot_agreement` must be a data.frame with a `spot` column.", call. = FALSE)
    }
    if (!("agreement" %in% colnames(dat))) {
      if ("mean_mad" %in% colnames(dat)) {
        dat$agreement <- 1 / (1 + dat$mean_mad)
      } else {
        stop("Spot agreement data must include `agreement` or `mean_mad`.", call. = FALSE)
      }
    }
    return(plot_spatial_metric_df(
      seu = x$seu,
      metrics = dat[, c("spot", "agreement"), drop = FALSE],
      value_col = "agreement",
      palette = palette,
      base_size = base_size,
      limits = c(0, 1),
      title = "Spot-Level Method Agreement",
      subtitle = "Higher values indicate stronger cross-method consistency",
      legend_title = "Agreement"
    ))
  }

  if (type %in% c("consensus_map", "confidence_map")) {
    return(plot_consensus_confidence(x, palette = palette, base_size = base_size))
  }
  if (type == "disagreement_map") {
    return(plot_disagreement_map(x, palette = palette, base_size = base_size))
  }
  if (type == "ranking") {
    return(plot_method_ranking(x, palette = palette, base_size = base_size))
  }

  stop("Unsupported `type`.", call. = FALSE)
}
