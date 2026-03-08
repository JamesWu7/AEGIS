#' Plot method comparison and consensus outputs
#'
#' @param x An `aegis` object.
#' @param type One of `heatmap`, `spot_agreement`, `consensus_map`.
#'
#' @return A ggplot object.
#' @export
plot_compare <- function(x, type = c("heatmap", "spot_agreement", "consensus_map")) {
  assert_is_aegis(x)
  type <- match.arg(type)

  if (type == "heatmap") {
    if (is.null(x$consensus$comparison$heatmap_matrix)) {
      stop("Comparison results not found. Run compare_methods() first.", call. = FALSE)
    }

    h <- as.data.frame(as.table(x$consensus$comparison$heatmap_matrix), stringsAsFactors = FALSE)
    colnames(h) <- c("method_1", "method_2", "agreement")

    p <- ggplot2::ggplot(h, ggplot2::aes(x = .data$method_1, y = .data$method_2, fill = .data$agreement)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0.5) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "Method agreement heatmap", x = NULL, y = NULL, fill = "Correlation")

    return(patchwork::wrap_plots(p))
  }

  coords <- get_spatial_coordinates(x$seu)

  if (type == "spot_agreement") {
    if (is.null(x$consensus$comparison$spot_agreement)) {
      stop("Comparison results not found. Run compare_methods() first.", call. = FALSE)
    }

    dat <- dplyr::left_join(x$consensus$comparison$spot_agreement, coords, by = c("spot" = "spot"))

    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$imagecol, y = .data$imagerow, color = .data$agreement)) +
      ggplot2::geom_point(size = 0.8) +
      ggplot2::scale_y_reverse() +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "Spot-level agreement", x = "imagecol", y = "imagerow", color = "Agreement")

    return(p)
  }

  if (type == "consensus_map") {
    if (is.null(x$consensus$result$spot_confidence)) {
      stop("Consensus result not found. Run compute_consensus() first.", call. = FALSE)
    }

    dat <- dplyr::left_join(x$consensus$result$spot_confidence, coords, by = c("spot" = "spot"))

    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$imagecol, y = .data$imagerow, color = .data$confidence)) +
      ggplot2::geom_point(size = 0.8) +
      ggplot2::scale_y_reverse() +
      ggplot2::coord_equal() +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = "Consensus confidence", x = "imagecol", y = "imagerow", color = "Confidence")

    return(p)
  }

  stop("Unsupported `type`.", call. = FALSE)
}
