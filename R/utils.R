# Internal utility helpers

#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
assert_is_seurat <- function(seu, arg = "seu") {
  if (!inherits(seu, "Seurat")) {
    stop(sprintf("`%s` must be a Seurat object.", arg), call. = FALSE)
  }
}

#' @keywords internal
assert_is_aegis <- function(x, arg = "x") {
  if (!inherits(x, "aegis")) {
    stop(sprintf("`%s` must be an aegis object. Use as_aegis() first.", arg), call. = FALSE)
  }
}

#' @keywords internal
coerce_to_numeric_matrix <- function(obj, method_name) {
  if (is.data.frame(obj)) {
    obj <- as.matrix(obj)
  } else if (!is.matrix(obj)) {
    obj <- tryCatch(
      as.matrix(obj),
      error = function(e) {
        stop(
          sprintf("Method '%s': value must be a matrix/data.frame or coercible to matrix.", method_name),
          call. = FALSE
        )
      }
    )
  }

  if (!is.matrix(obj)) {
    stop(
      sprintf("Method '%s': value must be coercible to matrix/data.frame.", method_name),
      call. = FALSE
    )
  }

  if (is.null(rownames(obj))) {
    stop(sprintf("Method '%s': row names are required.", method_name), call. = FALSE)
  }

  if (anyNA(rownames(obj)) || any(trimws(rownames(obj)) == "")) {
    stop(sprintf("Method '%s': missing/empty row names are not allowed.", method_name), call. = FALSE)
  }

  if (anyDuplicated(rownames(obj))) {
    stop(sprintf("Method '%s': duplicated row names are not allowed.", method_name), call. = FALSE)
  }
  if (nrow(obj) == 0L || ncol(obj) == 0L) {
    stop(sprintf("Method '%s': matrices with zero rows or zero columns are not allowed.", method_name), call. = FALSE)
  }
  if (is.null(colnames(obj))) {
    stop(sprintf("Method '%s': column names are required.", method_name), call. = FALSE)
  }
  if (anyNA(colnames(obj)) || any(trimws(colnames(obj)) == "")) {
    stop(sprintf("Method '%s': missing/empty column names are not allowed.", method_name), call. = FALSE)
  }
  if (anyDuplicated(colnames(obj))) {
    stop(sprintf("Method '%s': duplicated column names are not allowed.", method_name), call. = FALSE)
  }

  suppressWarnings(storage.mode(obj) <- "double")
  if (!is.numeric(obj)) {
    stop(sprintf("Method '%s': matrix must be numeric.", method_name), call. = FALSE)
  }
  if (all(is.na(obj))) {
    stop(sprintf("Method '%s': matrix coercion resulted in all NA values; numeric input is required.", method_name), call. = FALSE)
  }
  if (anyNA(obj)) {
    stop(sprintf("Method '%s': matrix contains NA values; numeric, non-missing values are required.", method_name), call. = FALSE)
  }

  obj
}

#' @keywords internal
align_matrix_to_spots <- function(mat, spots, method_name) {
  missing_spots <- setdiff(spots, rownames(mat))
  if (length(missing_spots) > 0) {
    preview <- paste(utils::head(missing_spots, 5L), collapse = ", ")
    stop(
      sprintf(
        "Method '%s': %d Seurat spots are missing from row names (e.g. %s).",
        method_name,
        length(missing_spots),
        preview
      ),
      call. = FALSE
    )
  }

  extra_spots <- setdiff(rownames(mat), spots)
  if (length(extra_spots) > 0) {
    preview <- paste(utils::head(extra_spots, 5L), collapse = ", ")
    stop(
      sprintf(
        "Method '%s': %d extra rows are not present in Seurat object (e.g. %s).",
        method_name,
        length(extra_spots),
        preview
      ),
      call. = FALSE
    )
  }

  if (!identical(rownames(mat), spots)) {
    # Explicitly reorder rows to match Seurat spot order.
    warning(
      sprintf("Method '%s': reordering rows to match Seurat spot order.", method_name),
      call. = FALSE
    )
    mat <- mat[spots, , drop = FALSE]
  }

  mat
}

#' @keywords internal
get_spatial_coordinates <- function(seu) {
  coords <- extract_spatial_coords(seu)
  out <- data.frame(
    spot = coords$spot,
    imagecol = coords$x,
    imagerow = coords$y,
    stringsAsFactors = FALSE
  )
  rownames(out) <- out$spot
  out
}

#' @keywords internal
extract_spatial_coords <- function(seu) {
  assert_is_seurat(seu)
  spots <- colnames(seu)
  if (is.null(spots) || length(spots) == 0L) {
    stop("Seurat object has no spot names (`colnames(seu)`).", call. = FALSE)
  }

  image_names <- names(seu@images)
  if (length(image_names) == 0L) {
    stop("No spatial image found in Seurat object (`seu@images` is empty).", call. = FALSE)
  }

  coords <- NULL
  for (img in image_names) {
    candidate <- tryCatch(
      Seurat::GetTissueCoordinates(seu, image = img),
      error = function(e) NULL
    )
    if (!is.null(candidate) && nrow(candidate) > 0L) {
      coords <- as.data.frame(candidate)
      break
    }
  }
  if (is.null(coords) || nrow(coords) == 0L) {
    stop("Could not extract spatial coordinates from Seurat image slots.", call. = FALSE)
  }

  if (!is.null(rownames(coords)) && all(nzchar(rownames(coords)))) {
    spot <- rownames(coords)
  } else if ("cell" %in% colnames(coords)) {
    spot <- as.character(coords$cell)
  } else if ("barcode" %in% colnames(coords)) {
    spot <- as.character(coords$barcode)
  } else {
    stop("Could not determine spot identifiers from tissue coordinates.", call. = FALSE)
  }

  x <- NULL
  y <- NULL
  if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
    x <- coords$imagecol
    y <- coords$imagerow
  } else if (all(c("x", "y") %in% colnames(coords))) {
    x <- coords$x
    y <- coords$y
  } else if (all(c("col", "row") %in% colnames(coords))) {
    x <- coords$col
    y <- coords$row
  } else if (all(c("pxl_col_in_fullres", "pxl_row_in_fullres") %in% colnames(coords))) {
    x <- coords$pxl_col_in_fullres
    y <- coords$pxl_row_in_fullres
  }
  if (is.null(x) || is.null(y)) {
    stop(
      "Spatial coordinates must include imagecol/imagerow, x/y, col/row, or pxl_col_in_fullres/pxl_row_in_fullres.",
      call. = FALSE
    )
  }

  out <- data.frame(
    spot = as.character(spot),
    x = as.numeric(x),
    y = as.numeric(y),
    stringsAsFactors = FALSE
  )
  if (anyNA(out$spot) || any(trimws(out$spot) == "")) {
    stop("Spatial coordinates include missing or empty spot IDs.", call. = FALSE)
  }
  if (anyDuplicated(out$spot)) {
    stop("Spatial coordinates include duplicated spot IDs.", call. = FALSE)
  }

  missing_spots <- setdiff(spots, out$spot)
  if (length(missing_spots) > 0L) {
    preview <- paste(utils::head(missing_spots, 5L), collapse = ", ")
    stop(sprintf("Could not match coordinates to all Seurat spots (e.g. %s).", preview), call. = FALSE)
  }

  out <- out[match(spots, out$spot), , drop = FALSE]
  rownames(out) <- out$spot
  out
}

#' @keywords internal
assemble_spatial_plot_data <- function(seu, metrics, spot_col = "spot") {
  if (!is.data.frame(metrics)) {
    stop("`metrics` must be a data.frame.", call. = FALSE)
  }
  if (!(spot_col %in% colnames(metrics))) {
    stop(sprintf("`metrics` must contain a `%s` column.", spot_col), call. = FALSE)
  }

  coords <- extract_spatial_coords(seu)
  dat <- dplyr::left_join(
    metrics,
    coords,
    by = stats::setNames("spot", spot_col)
  )
  if (anyNA(dat$x) || anyNA(dat$y)) {
    stop("Some spots in metrics could not be matched to spatial coordinates.", call. = FALSE)
  }
  dat
}

#' @keywords internal
get_plot_palette <- function(palette = "nature", type = c("continuous", "categorical", "diverging"), n = 256L) {
  type <- match.arg(type)
  palette <- tolower(palette)
  n <- max(3L, as.integer(n))

  nature_cont <- c("#0C2C84", "#225EA8", "#1D91C0", "#41B6C4", "#7FCDBB", "#C7E9B4", "#FFFFD9")
  nature_cat <- c("#2A6F97", "#D1495B", "#2A9D8F", "#E9C46A", "#264653", "#F4A261", "#8AB17D", "#7B6D8D")
  nature_div <- c("#3B4CC0", "#7B9FF9", "#C9D7F0", "#F7F7F7", "#F2B6A0", "#D6604D", "#B40426")

  viridis_cont <- grDevices::hcl.colors(n, "viridis")
  viridis_cat <- grDevices::hcl.colors(max(n, 8L), "viridis")[seq_len(max(8L, min(n, 12L)))]
  viridis_div <- grDevices::hcl.colors(n, "Purple-Green")

  scico_cont <- grDevices::hcl.colors(n, "Blue-Yellow 2")
  scico_cat <- grDevices::hcl.colors(max(n, 8L), "Dark 3")[seq_len(max(8L, min(n, 12L)))]
  scico_div <- grDevices::hcl.colors(n, "Blue-Red 3")

  brewer_cont <- grDevices::hcl.colors(n, "YlOrRd")
  brewer_cat <- grDevices::hcl.colors(max(n, 8L), "Set 2")[seq_len(max(8L, min(n, 12L)))]
  brewer_div <- grDevices::hcl.colors(n, "RdBu")

  values <- switch(palette,
    nature = switch(type, continuous = nature_cont, categorical = nature_cat, diverging = nature_div),
    viridis = switch(type, continuous = viridis_cont, categorical = viridis_cat, diverging = viridis_div),
    scico = switch(type, continuous = scico_cont, categorical = scico_cat, diverging = scico_div),
    brewer = switch(type, continuous = brewer_cont, categorical = brewer_cat, diverging = brewer_div),
    stop("Unsupported palette. Choose one of: nature, viridis, scico, brewer.", call. = FALSE)
  )

  if (type == "continuous" || type == "diverging") {
    grDevices::colorRampPalette(values)(n)
  } else {
    values
  }
}

#' @keywords internal
scale_fill_aegis <- function(palette = "nature", type = c("continuous", "categorical", "diverging"), na.value = "grey90", ...) {
  type <- match.arg(type)
  if (type == "categorical") {
    return(ggplot2::scale_fill_manual(values = get_plot_palette(palette, "categorical"), na.value = na.value, ...))
  }
  ggplot2::scale_fill_gradientn(colors = get_plot_palette(palette, type, n = 256L), na.value = na.value, ...)
}

#' @keywords internal
scale_color_aegis <- function(palette = "nature", type = c("continuous", "categorical", "diverging"), na.value = "grey90", ...) {
  type <- match.arg(type)
  if (type == "categorical") {
    return(ggplot2::scale_color_manual(values = get_plot_palette(palette, "categorical"), na.value = na.value, ...))
  }
  ggplot2::scale_color_gradientn(colors = get_plot_palette(palette, type, n = 256L), na.value = na.value, ...)
}

#' @keywords internal
theme_aegis <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size * 1.2, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.95),
      axis.title = ggplot2::element_text(size = base_size * 0.95),
      axis.text = ggplot2::element_text(size = base_size * 0.85, colour = "#222222"),
      legend.title = ggplot2::element_text(size = base_size * 0.9, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size * 0.8),
      strip.text = ggplot2::element_text(size = base_size * 0.9, face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "#D9D9D9"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key.height = grid::unit(0.9, "lines")
    )
}

#' @keywords internal
theme_aegis_spatial <- function(base_size = 12) {
  theme_aegis(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

#' @keywords internal
build_agreement_heatmap_data <- function(comparison) {
  dat <- comparison$pairwise_celltype_cor %||% comparison$pairwise_celltype
  if (is.null(dat) || !is.data.frame(dat) || nrow(dat) == 0L) {
    stop("Pairwise comparison data is unavailable for heatmap plotting.", call. = FALSE)
  }
  if (!all(c("method_1", "method_2", "celltype") %in% colnames(dat))) {
    stop("Pairwise comparison data must include method_1, method_2, and celltype columns.", call. = FALSE)
  }
  cor_col <- if ("correlation" %in% colnames(dat)) "correlation" else if ("pearson_cor" %in% colnames(dat)) "pearson_cor" else NULL
  if (is.null(cor_col)) {
    stop("Pairwise comparison data must include `correlation` or `pearson_cor`.", call. = FALSE)
  }
  dat <- dat |>
    dplyr::mutate(
      method_pair = paste(.data$method_1, .data$method_2, sep = " vs "),
      correlation = .data[[cor_col]]
    )
  dat
}

#' @keywords internal
format_report_table <- function(tbl, digits = 3, max_rows = 20L) {
  if (is.null(tbl) || !is.data.frame(tbl)) {
    return(data.frame())
  }
  out <- tbl
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], function(v) round(v, digits = digits))
  if (nrow(out) > max_rows) {
    out <- utils::head(out, max_rows)
  }
  out
}

#' @keywords internal
method_list <- function(x) {
  names(x$deconv)
}

#' @keywords internal
shared_celltypes_all_methods <- function(deconv_list) {
  Reduce(intersect, lapply(deconv_list, colnames))
}

#' @keywords internal
row_entropy <- function(mat, eps = 1e-12) {
  p <- pmax(mat, eps)
  -rowSums(p * log(p))
}

#' @keywords internal
nearest_neighbors <- function(coords, k = 6L) {
  n <- nrow(coords)
  if (n <= 1L) {
    return(vector("list", n))
  }

  k <- max(1L, min(as.integer(k), n - 1L))
  dmat <- as.matrix(stats::dist(coords[, c("imagecol", "imagerow")]))
  diag(dmat) <- Inf

  lapply(seq_len(n), function(i) {
    order(dmat[i, ], decreasing = FALSE)[seq_len(k)]
  })
}

#' @keywords internal
resolve_template <- function() {
  pkg_template <- system.file("templates", "aegis_report.Rmd", package = "AEGIS")
  if (nzchar(pkg_template) && file.exists(pkg_template)) {
    return(pkg_template)
  }

  local_template <- file.path("inst", "templates", "aegis_report.Rmd")
  if (file.exists(local_template)) {
    return(normalizePath(local_template, winslash = "/", mustWork = TRUE))
  }

  stop("Could not locate report template `inst/templates/aegis_report.Rmd`.", call. = FALSE)
}

#' @keywords internal
resolve_path_in_dir <- function(data_dir, filename, required = TRUE) {
  path <- file.path(data_dir, filename)
  if (!file.exists(path) && required) {
    stop(sprintf("Required file not found: %s", normalizePath(path, winslash = "/", mustWork = FALSE)), call. = FALSE)
  }
  path
}

utils::globalVariables(c(".data"))
