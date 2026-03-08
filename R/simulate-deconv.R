#' Simulate deconvolution results for a Seurat spatial object
#'
#' Generates realistic mock spot-by-celltype proportion matrices for multiple
#' methods. This is intended for MVP development and demos before plugging in
#' real deconvolution backends.
#'
#' @param seu A Seurat object.
#' @param cell_types Character vector of cell types.
#' @param methods Character vector of method names.
#' @param seed Random seed for reproducibility.
#' @param noise_scale Numeric scale for method-specific perturbation.
#'
#' @return A named list of numeric matrices (spot-by-celltype proportions).
#' @export
simulate_deconv_results <- function(
    seu,
    cell_types = c("B_cell", "T_cell", "NK_cell", "Plasma", "Myeloid", "Endothelial", "Stromal"),
    methods = c("RCTD", "SPOTlight", "cell2location"),
    seed = 123,
    noise_scale = 0.05
) {
  assert_is_seurat(seu, "seu")

  spots <- colnames(seu)
  if (is.null(spots) || length(spots) == 0L) {
    stop("`seu` must contain spot/cell column names.", call. = FALSE)
  }
  if (anyNA(spots) || any(trimws(spots) == "") || anyDuplicated(spots)) {
    stop("`seu` spot names must be non-empty and unique.", call. = FALSE)
  }

  if (!is.character(cell_types) || length(cell_types) < 2L || any(trimws(cell_types) == "")) {
    stop("`cell_types` must be a character vector with at least 2 non-empty entries.", call. = FALSE)
  }
  if (anyDuplicated(cell_types)) {
    stop("`cell_types` must be unique.", call. = FALSE)
  }
  if (!is.character(methods) || length(methods) == 0L || any(trimws(methods) == "")) {
    stop("`methods` must be a non-empty character vector of method names.", call. = FALSE)
  }
  if (anyDuplicated(methods)) {
    stop("`methods` must be unique.", call. = FALSE)
  }
  if (!is.numeric(noise_scale) || length(noise_scale) != 1L || is.na(noise_scale) || noise_scale < 0) {
    stop("`noise_scale` must be a single non-negative numeric value.", call. = FALSE)
  }

  n_spots <- length(spots)
  n_ct <- length(cell_types)

  set.seed(seed)
  coords <- tryCatch(get_spatial_coordinates(seu), error = function(e) NULL)
  if (!is.null(coords)) {
    coords <- coords[spots, , drop = FALSE]
    grad_x <- as.numeric(scale(coords$imagecol))
    grad_y <- as.numeric(scale(coords$imagerow))
  } else {
    grad_x <- as.numeric(scale(stats::rnorm(n_spots)))
    grad_y <- as.numeric(scale(stats::rnorm(n_spots)))
  }

  # Build a latent composition with spatial gradients + spot dominance.
  base_alpha <- seq(1.1, 2.5, length.out = n_ct)
  latent <- matrix(0, nrow = n_spots, ncol = n_ct, dimnames = list(spots, cell_types))
  for (j in seq_len(n_ct)) {
    phase <- (j - 1) * pi / (n_ct + 1)
    spatial_signal <- 0.45 * cos(phase) * grad_x + 0.35 * sin(phase) * grad_y
    alpha_j <- pmax(0.12, base_alpha[[j]] * exp(spatial_signal))
    latent[, j] <- stats::rgamma(n_spots, shape = alpha_j, rate = 1)
  }
  latent <- latent / rowSums(latent)

  dominant_frac <- 0.35
  n_dom <- max(1L, floor(n_spots * dominant_frac))
  dominant_idx <- sample.int(n_spots, size = n_dom)
  dominant_ct <- sample.int(n_ct, size = n_dom, replace = TRUE)
  boost <- stats::runif(n_dom, min = 1.6, max = 3.0)
  latent[cbind(dominant_idx, dominant_ct)] <- latent[cbind(dominant_idx, dominant_ct)] * boost
  latent <- latent / rowSums(latent)

  center_profile <- matrix(colMeans(latent), nrow = n_spots, ncol = n_ct, byrow = TRUE)

  method_profiles <- list(
    RCTD = list(power = 1.20, center_mix = 0.03, noise_mult = 0.80),
    SPOTlight = list(power = 0.88, center_mix = 0.12, noise_mult = 1.00),
    cell2location = list(power = 1.02, center_mix = 0.08, noise_mult = 1.15)
  )

  default_profile <- list(power = 1.00, center_mix = 0.06, noise_mult = 1.00)

  make_method_matrix <- function(method) {
    profile <- method_profiles[[method]]
    if (is.null(profile)) {
      profile <- default_profile
    }

    mat <- latent ^ profile$power
    mat <- (1 - profile$center_mix) * mat + profile$center_mix * center_profile
    mat <- mat + matrix(
      stats::rnorm(n_spots * n_ct, mean = 0, sd = noise_scale * profile$noise_mult),
      nrow = n_spots,
      ncol = n_ct
    )
    mat <- pmax(mat, 0)
    rs <- rowSums(mat)
    rs[rs <= 0] <- 1
    mat <- mat / rs
    rownames(mat) <- spots
    colnames(mat) <- cell_types
    mat
  }

  out <- lapply(methods, make_method_matrix)
  names(out) <- methods

  # Post-simulation validation for deterministic and inspectable outputs.
  for (method in names(out)) {
    mat <- out[[method]]
    if (!is.matrix(mat) || !is.numeric(mat)) {
      stop(sprintf("Simulation output for '%s' is not a numeric matrix.", method), call. = FALSE)
    }
    if (anyNA(mat)) {
      stop(sprintf("Simulation output for '%s' contains NA values.", method), call. = FALSE)
    }
    if (any(mat < 0)) {
      stop(sprintf("Simulation output for '%s' contains negative values.", method), call. = FALSE)
    }
    if (!identical(rownames(mat), spots)) {
      stop(sprintf("Simulation output for '%s' rownames do not match Seurat spots.", method), call. = FALSE)
    }
    if (!identical(colnames(mat), cell_types)) {
      stop(sprintf("Simulation output for '%s' colnames do not match `cell_types`.", method), call. = FALSE)
    }
    if (any(abs(rowSums(mat) - 1) > 1e-6)) {
      stop(sprintf("Simulation output for '%s' has row sums that deviate from 1.", method), call. = FALSE)
    }
  }

  out
}
