#' Simulate deconvolution results for a Seurat spatial object
#'
#' Generates realistic mock spot-by-celltype proportion matrices for multiple
#' methods. This is intended for MVP development and demos before plugging in
#' real deconvolution backends.
#'
#' @param seu A Seurat object.
#' @param methods Character vector of method names.
#' @param cell_types Character vector of cell types.
#' @param seed Random seed for reproducibility.
#'
#' @return A named list of numeric matrices (spot-by-celltype proportions).
#' @export
simulate_deconv_results <- function(
    seu,
    methods = c("RCTD", "SPOTlight", "cell2location"),
    cell_types = c("B_cell", "T_cell", "NK_cell", "Plasma", "Myeloid", "Endothelial", "Stromal"),
    seed = 123
) {
  assert_is_seurat(seu, "seu")

  if (length(methods) < 2L) {
    stop("`methods` must contain at least two method names.", call. = FALSE)
  }
  if (anyDuplicated(methods) || any(trimws(methods) == "")) {
    stop("`methods` must be unique and non-empty.", call. = FALSE)
  }

  spots <- colnames(seu)
  n_spots <- length(spots)
  n_ct <- length(cell_types)
  if (n_spots == 0L || n_ct < 2L) {
    stop("Seurat object must contain spots and `cell_types` must have at least 2 entries.", call. = FALSE)
  }

  set.seed(seed)

  coords <- tryCatch(get_spatial_coordinates(seu), error = function(e) NULL)
  if (!is.null(coords)) {
    coords <- coords[spots, , drop = FALSE]
    x <- as.numeric(scale(coords$imagecol))
    y <- as.numeric(scale(coords$imagerow))
  } else {
    x <- stats::rnorm(n_spots)
    y <- stats::rnorm(n_spots)
  }

  # Base abundance tendencies across cell types.
  base_alpha <- seq(1.6, 3.2, length.out = n_ct)

  # Spatial effects create plausible local gradients for selected populations.
  spatial_effect <- matrix(0, nrow = n_spots, ncol = n_ct)
  colnames(spatial_effect) <- cell_types

  idx <- function(ct) which(cell_types == ct)
  if (length(idx("T_cell")) > 0L) spatial_effect[, idx("T_cell")] <- 0.35 * x
  if (length(idx("B_cell")) > 0L) spatial_effect[, idx("B_cell")] <- -0.25 * x
  if (length(idx("Myeloid")) > 0L) spatial_effect[, idx("Myeloid")] <- 0.30 * y
  if (length(idx("Endothelial")) > 0L) spatial_effect[, idx("Endothelial")] <- -0.28 * y
  if (length(idx("Stromal")) > 0L) spatial_effect[, idx("Stromal")] <- 0.20 * (x + y)

  latent <- matrix(0, nrow = n_spots, ncol = n_ct)
  for (j in seq_len(n_ct)) {
    alpha_j <- pmax(0.2, base_alpha[[j]] * exp(spatial_effect[, j]))
    latent[, j] <- stats::rgamma(n_spots, shape = alpha_j, rate = 1)
  }
  latent <- latent / rowSums(latent)

  # Method-specific perturbation profiles.
  method_bias <- lapply(methods, function(m) rep(1, n_ct))
  names(method_bias) <- methods

  if ("RCTD" %in% methods && n_ct >= 3L) {
    b <- rep(1, n_ct)
    b[cell_types %in% c("T_cell", "B_cell")] <- 1.08
    b[cell_types %in% c("Stromal")] <- 0.94
    method_bias[["RCTD"]] <- b
  }
  if ("SPOTlight" %in% methods && n_ct >= 3L) {
    b <- rep(1, n_ct)
    b[cell_types %in% c("Myeloid", "Plasma")] <- 1.10
    b[cell_types %in% c("NK_cell")] <- 0.92
    method_bias[["SPOTlight"]] <- b
  }
  if ("cell2location" %in% methods && n_ct >= 3L) {
    b <- rep(1, n_ct)
    b[cell_types %in% c("Endothelial", "Stromal")] <- 1.12
    b[cell_types %in% c("Plasma")] <- 0.90
    method_bias[["cell2location"]] <- b
  }

  out <- lapply(methods, function(method) {
    noise_sd <- switch(method,
      "RCTD" = 0.03,
      "SPOTlight" = 0.05,
      "cell2location" = 0.04,
      0.04
    )

    mat <- latent * matrix(method_bias[[method]], nrow = n_spots, ncol = n_ct, byrow = TRUE)
    mat <- mat + matrix(stats::rnorm(n_spots * n_ct, mean = 0, sd = noise_sd), n_spots, n_ct)
    mat <- pmax(mat, 0)

    rs <- rowSums(mat)
    rs[rs <= 0] <- 1
    mat <- mat / rs

    rownames(mat) <- spots
    colnames(mat) <- cell_types
    mat
  })

  names(out) <- methods
  out
}
