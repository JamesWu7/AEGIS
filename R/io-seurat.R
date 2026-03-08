#' Load the 10x Human Lymph Node Visium dataset as a Seurat object
#'
#' This loader is tailored to the authoritative MVP assets:
#' `V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5`,
#' `V1_Human_Lymph_Node_spatial.tar.gz`, and
#' `V1_Human_Lymph_Node_metrics_summary.csv`.
#'
#' It will unpack the spatial tarball automatically when `spatial/` is missing,
#' then call [Seurat::Load10X_Spatial()].
#'
#' @param data_dir Directory containing the raw files.
#' @param h5_file Name of the feature-barcode H5 file.
#' @param spatial_tar Name of the spatial tar.gz archive.
#' @param metrics_file Name of the metrics CSV file.
#' @param assay Assay name passed to [Seurat::Load10X_Spatial()].
#' @param slice Slice name passed to [Seurat::Load10X_Spatial()].
#' @param filter.matrix Logical; forwarded to [Seurat::Load10X_Spatial()].
#'
#' @return A Seurat spatial object.
#' @export
load_10x_lymphnode <- function(
    data_dir = ".",
    h5_file = "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5",
    spatial_tar = "V1_Human_Lymph_Node_spatial.tar.gz",
    metrics_file = "V1_Human_Lymph_Node_metrics_summary.csv",
    assay = "Spatial",
    slice = "slice1",
    filter.matrix = TRUE
) {
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(data_dir)) {
    stop(sprintf("`data_dir` does not exist: %s", data_dir), call. = FALSE)
  }

  h5_path <- resolve_path_in_dir(data_dir, h5_file, required = TRUE)
  metrics_path <- resolve_path_in_dir(data_dir, metrics_file, required = FALSE)
  if (!file.exists(metrics_path)) {
    warning(
      sprintf(
        "Metrics file not found at %s. Continuing because metrics are not required by Load10X_Spatial.",
        normalizePath(metrics_path, winslash = "/", mustWork = FALSE)
      ),
      call. = FALSE
    )
  }

  spatial_dir <- file.path(data_dir, "spatial")
  if (!is_valid_spatial_dir(spatial_dir)) {
    tar_path <- resolve_path_in_dir(data_dir, spatial_tar, required = TRUE)
    spatial_dir <- ensure_spatial_dir(data_dir = data_dir, tar_path = tar_path)
  }

  if (!is_valid_spatial_dir(spatial_dir)) {
    stop(
      sprintf(
        "Malformed spatial directory at %s. Expected `scalefactors_json.json` and `tissue_positions*.csv`.",
        spatial_dir
      ),
      call. = FALSE
    )
  }

  seu <- Seurat::Load10X_Spatial(
    data.dir = data_dir,
    filename = basename(h5_path),
    assay = assay,
    slice = slice,
    filter.matrix = filter.matrix
  )

  seu
}

#' @keywords internal
is_valid_spatial_dir <- function(spatial_dir) {
  if (!dir.exists(spatial_dir)) {
    return(FALSE)
  }

  scale_ok <- file.exists(file.path(spatial_dir, "scalefactors_json.json"))
  tissue_ok <- length(list.files(spatial_dir, pattern = "^tissue_positions(_list)?\\.csv$")) > 0L
  image_ok <- any(file.exists(file.path(spatial_dir, c("tissue_lowres_image.png", "tissue_hires_image.png"))))

  scale_ok && tissue_ok && image_ok
}

#' @keywords internal
ensure_spatial_dir <- function(data_dir, tar_path) {
  unpack_dir <- tempfile("aegis_spatial_unpack_")
  dir.create(unpack_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(unpack_dir, recursive = TRUE, force = TRUE), add = TRUE)

  utils::untar(tar_path, exdir = unpack_dir)

  candidate <- find_spatial_candidate(unpack_dir)
  if (is.null(candidate)) {
    stop(
      sprintf("Could not find a valid `spatial/` directory after unpacking: %s", tar_path),
      call. = FALSE
    )
  }

  target <- file.path(data_dir, "spatial")
  if (dir.exists(target)) {
    unlink(target, recursive = TRUE, force = TRUE)
  }
  dir.create(target, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(candidate, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  copied <- file.copy(files, target, recursive = TRUE)
  if (!all(copied)) {
    stop("Failed to materialize `spatial/` directory from tarball.", call. = FALSE)
  }

  target
}

#' @keywords internal
find_spatial_candidate <- function(root_dir) {
  direct <- file.path(root_dir, "spatial")
  if (is_valid_spatial_dir(direct)) {
    return(direct)
  }

  dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
  for (d in dirs) {
    if (is_valid_spatial_dir(d)) {
      return(d)
    }
  }

  NULL
}
