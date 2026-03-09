#' Run spatial smoothness audit
#'
#' Uses nearest-neighbor neighborhoods on tissue coordinates to estimate local
#' inconsistency of predicted cell-type proportions.
#'
#' @param x An `aegis` object.
#' @param k Number of nearest neighbors.
#'
#' @return The modified `aegis` object with `x$audit$spatial` populated.
#' @export
audit_spatial <- function(x, k = 6L) {
  if (is_multi_sample_context(x)) {
    sample_results <- iterate_aegis_samples(x, function(obj) audit_spatial(obj, k = k))
    by_sample <- lapply(sample_results, function(obj) obj$audit$spatial)

    spot_metrics <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]$spot_metrics %||% data.frame()
      if (nrow(tbl) == 0L) return(tbl)
      tbl$sample_id <- sid
      tbl
    }))
    detail <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]$detail %||% data.frame()
      if (nrow(tbl) == 0L) return(tbl)
      tbl$sample_id <- sid
      tbl
    }))
    method_summary <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]$method_summary %||% data.frame()
      if (nrow(tbl) == 0L) return(tbl)
      tbl$sample_id <- sid
      tbl
    }))
    summary <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]$summary %||% data.frame()
      if (nrow(tbl) == 0L) return(tbl)
      tbl$sample_id <- sid
      tbl
    }))

    x$audit$spatial <- list(
      by_sample = by_sample,
      spot_metrics = spot_metrics,
      detail = detail,
      method_summary = method_summary,
      summary = summary,
      params = list(k = k)
    )
    return(x)
  }

  assert_is_aegis(x)
  assert_is_seurat(x$seu, "x$seu")
  if (is.null(x$deconv) || !is.list(x$deconv) || length(x$deconv) == 0L) {
    stop("`x$deconv` must be a non-empty list.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0 || k %% 1 != 0) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }
  k <- as.integer(k)

  coords <- get_spatial_coordinates(x$seu)
  spots <- colnames(x$seu)
  if (!all(spots %in% coords$spot)) {
    missing_spots <- setdiff(spots, coords$spot)
    preview <- paste(utils::head(missing_spots, 5L), collapse = ", ")
    stop(sprintf("Could not match spatial coordinates to all spots (e.g. %s).", preview), call. = FALSE)
  }
  coords <- coords[match(spots, coords$spot), , drop = FALSE]
  rownames(coords) <- spots

  n_spots <- length(spots)
  if (n_spots < 2L) {
    stop("At least two spots are required for spatial audit.", call. = FALSE)
  }
  if (k >= n_spots) {
    stop(sprintf("`k` must be <= n_spots - 1 (%d).", n_spots - 1L), call. = FALSE)
  }

  nn <- nearest_neighbors(coords, k = k)

  spot_tables <- lapply(names(x$deconv), function(method) {
    mat <- x$deconv[[method]]
    if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
      stop(sprintf("Method '%s' must be a non-empty numeric matrix.", method), call. = FALSE)
    }
    if (!identical(rownames(mat), spots)) {
      stop(sprintf("Method '%s' rownames are not aligned to Seurat spots.", method), call. = FALSE)
    }
    ct_names <- colnames(mat)

    ct_tables <- lapply(seq_along(ct_names), function(j) {
      v <- mat[, j]
      local_dev <- vapply(seq_along(v), function(i) {
        nbr <- nn[[i]]
        abs(v[[i]] - mean(v[nbr]))
      }, numeric(1))

      data.frame(
        spot = spots,
        method = method,
        celltype = ct_names[[j]],
        local_inconsistency = local_dev,
        stringsAsFactors = FALSE
      )
    })

    dplyr::bind_rows(ct_tables)
  })

  detail <- dplyr::bind_rows(spot_tables)

  spot_metrics <- detail |>
    dplyr::group_by(.data$spot, .data$method) |>
    dplyr::summarise(
      local_inconsistency = mean(.data$local_inconsistency, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(smoothness = 1 / (1 + .data$local_inconsistency))

  method_celltype_summary <- detail |>
    dplyr::group_by(.data$method, .data$celltype) |>
    dplyr::summarise(
      mean_local_inconsistency = mean(.data$local_inconsistency, na.rm = TRUE),
      .groups = "drop"
    )

  summary <- spot_metrics |>
    dplyr::group_by(.data$method) |>
    dplyr::summarise(
      mean_local_inconsistency = mean(.data$local_inconsistency, na.rm = TRUE),
      mean_smoothness = mean(.data$smoothness, na.rm = TRUE),
      .groups = "drop"
    )

  names(nn) <- spots

  x$audit$spatial <- list(
    spot_metrics = spot_metrics,
    detail = detail,
    spot_summary = spot_metrics,
    method_summary = method_celltype_summary,
    summary = summary,
    neighbors = nn,
    params = list(k = k)
  )

  x
}
