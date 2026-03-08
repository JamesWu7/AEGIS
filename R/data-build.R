#' Default marker list used by AEGIS MVP
#'
#' @return A named list of marker genes.
#' @export
aegis_default_markers <- function() {
  list(
    B_cell = c("CD79A", "MS4A1", "CD74"),
    T_cell = c("CD3D", "CD3E", "TRBC1"),
    NK_cell = c("NKG7", "GNLY", "KLRD1"),
    Plasma = c("MZB1", "JCHAIN", "SDC1"),
    Myeloid = c("LYZ", "FCER1G", "TYMP"),
    Endothelial = c("PECAM1", "VWF", "KDR"),
    Stromal = c("COL1A1", "COL1A2", "DCN")
  )
}

#' Build and save example data from raw Human Lymph Node files
#'
#' @param raw_dir Directory containing authoritative raw files.
#' @param output_file Path to `.rda` output file.
#' @param object_name Name of the object saved in `.rda`.
#' @param max_spots Maximum number of spots to keep for package-size control.
#' @param seed Random seed for reproducible spot subsetting.
#'
#' @return Invisibly returns the saved Seurat object.
#' @export
build_aegis_example_data <- function(
    raw_dir = ".",
    output_file = file.path("data", "aegis_example.rda"),
    object_name = "aegis_example",
    max_spots = 1200L,
    seed = 42
) {
  seu <- load_10x_lymphnode(data_dir = raw_dir)

  set.seed(seed)
  all_spots <- colnames(seu)
  if (length(all_spots) > max_spots) {
    keep <- sort(sample(all_spots, size = max_spots))
    seu <- subset(seu, cells = keep)
  }

  if ("Spatial" %in% names(seu@assays)) {
    seu <- Seurat::NormalizeData(seu, assay = "Spatial", verbose = FALSE)
  }

  # Keep the object compact while preserving spatial plotting and marker scoring.
  seu <- Seurat::DietSeurat(
    seu,
    assays = "Spatial",
    dimreducs = NULL,
    graphs = NULL,
    misc = TRUE
  )

  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  assign(object_name, seu)
  save(list = object_name, file = output_file, compress = "bzip2")

  invisible(seu)
}
