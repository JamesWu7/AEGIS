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
  assert_is_aegis(x)

  coords <- get_spatial_coordinates(x$seu)
  coords <- coords[colnames(x$seu), , drop = FALSE]

  nn <- nearest_neighbors(coords, k = k)
  spots <- colnames(x$seu)

  spot_tables <- lapply(names(x$deconv), function(method) {
    mat <- x$deconv[[method]]
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

  spot_summary <- detail |>
    dplyr::group_by(.data$spot, .data$method) |>
    dplyr::summarise(
      smoothness_score = 1 / (1 + mean(.data$local_inconsistency, na.rm = TRUE)),
      mean_local_inconsistency = mean(.data$local_inconsistency, na.rm = TRUE),
      .groups = "drop"
    )

  method_summary <- detail |>
    dplyr::group_by(.data$method, .data$celltype) |>
    dplyr::summarise(
      mean_local_inconsistency = mean(.data$local_inconsistency, na.rm = TRUE),
      .groups = "drop"
    )

  x$audit$spatial <- list(
    detail = detail,
    spot_summary = spot_summary,
    method_summary = method_summary,
    params = list(k = k)
  )

  x
}
