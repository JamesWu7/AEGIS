#' Plot audit outputs
#'
#' @param x An `aegis` object.
#' @param type One of `sumdev`, `dominance`, `entropy`, `marker`, `smoothness`.
#' @param method Optional method name to subset to one deconvolution method.
#' @param sample Optional sample ID for multi-sample objects.
#' @param cell_type Optional cell type for marker/smoothness views.
#' @param palette Palette family: `nature`, `viridis`, `scico`, or `brewer`.
#' @param point_size Optional spatial point size. Auto-selected when `NULL`.
#' @param base_size Base font size for the plot theme.
#'
#' @return A ggplot object.
#' @export
plot_audit <- function(
    x,
    type = c("sumdev", "dominance", "entropy", "marker", "smoothness"),
    method = NULL,
    sample = NULL,
    cell_type = NULL,
    palette = "nature",
    point_size = NULL,
    base_size = 12) {
  x <- resolve_plot_input_aegis(x, sample = sample, fun_name = "plot_audit")

  assert_is_aegis(x)
  type <- match.arg(type)

  if (is.null(point_size)) {
    point_size <- auto_spatial_point_size(ncol(x$seu))
  }

  if (type %in% c("sumdev", "dominance", "entropy")) {
    if (is.null(x$audit$basic$spot_metrics)) {
      stop("Basic audit not found. Run audit_basic() first.", call. = FALSE)
    }
    metric_col <- switch(type,
      sumdev = "sum_dev",
      dominance = "dominance",
      entropy = "entropy"
    )

    dat <- x$audit$basic$spot_metrics
    if (!is.null(method)) {
      methods <- unique(dat$method)
      if (!(method %in% methods)) {
        stop(sprintf("Method '%s' not found in basic audit results.", method), call. = FALSE)
      }
      dat <- dat[dat$method == method, , drop = FALSE]
    }
    plot_dat <- dat[, c("spot", "method", metric_col), drop = FALSE]
    return(plot_spatial_metric_df(
      seu = x$seu,
      metrics = plot_dat,
      value_col = metric_col,
      palette = palette,
      base_size = base_size,
      point_size = point_size,
      facet_col = if (is.null(method) && length(unique(plot_dat$method)) > 1L) "method" else NULL,
      title = switch(type,
        sumdev = "Row-Sum Deviation Map",
        dominance = "Dominance Map",
        entropy = "Entropy Map"
      ),
      subtitle = if (!is.null(method)) paste("Method:", method) else "Faceted by method",
      legend_title = switch(type,
        sumdev = "Sum deviation",
        dominance = "Dominance",
        entropy = "Entropy"
      )
    ))
  }

  if (type == "marker") {
    if (is.null(x$audit$marker$concordance)) {
      stop("Marker audit not found. Run audit_marker() first.", call. = FALSE)
    }
    dat <- x$audit$marker$concordance
    if (!is.data.frame(dat) || nrow(dat) == 0L) {
      stop("Marker concordance table is empty.", call. = FALSE)
    }

    y_col <- if ("pearson_cor" %in% colnames(dat)) {
      "pearson_cor"
    } else if ("correlation" %in% colnames(dat)) {
      "correlation"
    } else {
      stop("Marker concordance table must contain `pearson_cor` or `correlation`.", call. = FALSE)
    }

    if (!is.null(method)) {
      if (!(method %in% unique(dat$method))) {
        stop(sprintf("Method '%s' not found in marker audit results.", method), call. = FALSE)
      }
      dat <- dat[dat$method == method, , drop = FALSE]
    }
    if (!is.null(cell_type)) {
      if (!(cell_type %in% unique(dat$celltype))) {
        stop(sprintf("Cell type '%s' not found in marker audit results.", cell_type), call. = FALSE)
      }
      dat <- dat[dat$celltype == cell_type, , drop = FALSE]
    }
    dat <- dat[is.finite(dat[[y_col]]), , drop = FALSE]
    if (nrow(dat) == 0L) {
      stop("No finite marker concordance values available to plot.", call. = FALSE)
    }

    dat$celltype <- stats::reorder(dat$celltype, dat[[y_col]], FUN = mean, na.rm = TRUE)

    p <- ggplot2::ggplot(
      dat,
      ggplot2::aes(x = .data$celltype, y = .data$method, color = .data[[y_col]], size = .data$n_markers_used)
    ) +
      ggplot2::geom_point(alpha = 0.95) +
      scale_color_aegis(palette = palette, type = "diverging", limits = c(-1, 1)) +
      ggplot2::scale_size_continuous(range = c(2.2, 6), breaks = c(1, 2, 3), limits = c(0, NA)) +
      theme_aegis(base_size = base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)) +
      ggplot2::labs(
        title = "Marker Concordance Summary",
        subtitle = "Method-celltype correlation between marker score and predicted abundance",
        x = "Cell type",
        y = "Method",
        color = "Correlation",
        size = "Markers used"
      )

    return(p)
  }

  if (type == "smoothness") {
    if (is.null(x$audit$spatial$spot_metrics)) {
      stop("Spatial audit not found. Run audit_spatial() first.", call. = FALSE)
    }

    dat <- x$audit$spatial$spot_metrics
    if (!all(c("spot", "method") %in% colnames(dat))) {
      stop("Spatial spot metrics must include `spot` and `method` columns.", call. = FALSE)
    }

    if (!is.null(cell_type)) {
      detail <- x$audit$spatial$detail
      if (is.null(detail) || !is.data.frame(detail) || nrow(detail) == 0L) {
        stop("Spatial detail metrics are unavailable for cell-type smoothness plotting.", call. = FALSE)
      }
      if (!(cell_type %in% unique(detail$celltype))) {
        stop(sprintf("Cell type '%s' not found in spatial audit detail results.", cell_type), call. = FALSE)
      }
      dat <- detail[detail$celltype == cell_type, c("spot", "method", "local_inconsistency"), drop = FALSE]
      dat$smoothness <- 1 / (1 + dat$local_inconsistency)
    } else {
      if (!("smoothness" %in% colnames(dat))) {
        if (!("local_inconsistency" %in% colnames(dat))) {
          stop("Spatial spot metrics must contain `smoothness` or `local_inconsistency`.", call. = FALSE)
        }
        dat$smoothness <- 1 / (1 + dat$local_inconsistency)
      }
      dat <- dat[, c("spot", "method", "smoothness"), drop = FALSE]
    }

    if (!is.null(method)) {
      if (!(method %in% unique(dat$method))) {
        stop(sprintf("Method '%s' not found in spatial audit results.", method), call. = FALSE)
      }
      dat <- dat[dat$method == method, , drop = FALSE]
    }

    plot_dat <- dat[, c("spot", "method", "smoothness"), drop = FALSE]
    return(plot_spatial_metric_df(
      seu = x$seu,
      metrics = plot_dat,
      value_col = "smoothness",
      palette = palette,
      base_size = base_size,
      point_size = point_size,
      facet_col = if (is.null(method) && length(unique(plot_dat$method)) > 1L) "method" else NULL,
      title = "Spatial Smoothness Map",
      subtitle = if (!is.null(cell_type)) paste("Cell type:", cell_type) else "Spot-level smoothness",
      legend_title = "Smoothness"
    ))
  }

  stop("Unsupported `type`.", call. = FALSE)
}
