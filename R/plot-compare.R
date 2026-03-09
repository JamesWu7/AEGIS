#' Plot method comparison and consensus outputs
#'
#' @param x An `aegis` object.
#' @param type One of `heatmap`, `spot_agreement`, `consensus_map`.
#' @param sample Optional sample ID for multi-sample objects.
#' @param palette Palette family: `nature`, `viridis`, `scico`, or `brewer`.
#' @param base_size Base font size for the plot theme.
#'
#' @return A ggplot object.
#' @export
plot_compare <- function(
    x,
    type = c("heatmap", "spot_agreement", "consensus_map"),
    sample = NULL,
    palette = "nature",
    base_size = 12) {
  if (is_multi_sample_context(x)) {
    sample_objs <- split_aegis_by_sample(x)
    if (is.null(sample)) {
      stop("Multi-sample object detected. Please provide `sample` for `plot_compare()`.", call. = FALSE)
    }
    if (!(sample %in% names(sample_objs))) {
      stop(sprintf("Sample '%s' not found. Available: %s", sample, paste(names(sample_objs), collapse = ", ")), call. = FALSE)
    }
    return(plot_compare(
      x = sample_objs[[sample]],
      type = type,
      sample = NULL,
      palette = palette,
      base_size = base_size
    ))
  }

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
    dat$method_pair <- factor(dat$method_pair, levels = pair_order)
    dat$celltype <- factor(dat$celltype, levels = sort(unique(dat$celltype)))

    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$method_pair, y = .data$celltype, fill = .data$correlation)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_aegis(palette = palette, type = "diverging", limits = c(-1, 1)) +
      theme_aegis(base_size = base_size) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1)
      ) +
      ggplot2::labs(
        title = "Method-Pair Agreement by Cell Type",
        subtitle = "Correlation across shared spots",
        x = "Method pair",
        y = "Cell type",
        fill = "Correlation"
      )
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
    spatial_dat <- assemble_spatial_plot_data(x$seu, dat[, c("spot", "agreement"), drop = FALSE])

    p <- ggplot2::ggplot(spatial_dat, ggplot2::aes(x = .data$x, y = .data$y, color = .data$agreement)) +
      ggplot2::geom_point(size = if (nrow(spatial_dat) > 5000L) 0.45 else 0.7, alpha = 0.95) +
      ggplot2::coord_fixed() +
      ggplot2::scale_y_reverse() +
      scale_color_aegis(palette = palette, type = "continuous", limits = c(0, 1)) +
      theme_aegis_spatial(base_size = base_size) +
      ggplot2::labs(
        title = "Spot-Level Method Agreement",
        subtitle = "Higher values indicate stronger cross-method consistency",
        color = "Agreement"
      )
    return(p)
  }

  if (is.null(x$consensus$result$spot_confidence)) {
    stop("Consensus result not found. Run compute_consensus() first.", call. = FALSE)
  }
  dat <- x$consensus$result$spot_confidence
  if (!is.data.frame(dat) || !all(c("spot", "confidence") %in% colnames(dat))) {
    stop("`x$consensus$result$spot_confidence` must include `spot` and `confidence`.", call. = FALSE)
  }
  spatial_dat <- assemble_spatial_plot_data(x$seu, dat[, c("spot", "confidence"), drop = FALSE])

  p <- ggplot2::ggplot(spatial_dat, ggplot2::aes(x = .data$x, y = .data$y, color = .data$confidence)) +
    ggplot2::geom_point(size = if (nrow(spatial_dat) > 5000L) 0.45 else 0.7, alpha = 0.95) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    scale_color_aegis(palette = palette, type = "continuous", limits = c(0, 1)) +
    theme_aegis_spatial(base_size = base_size) +
    ggplot2::labs(
      title = "Consensus Confidence Map",
      subtitle = "Per-spot confidence derived from cross-method disagreement",
      color = "Confidence"
    )

  p
}
