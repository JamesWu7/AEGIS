#!/usr/bin/env Rscript

# Reproducible data build entrypoint for AEGIS MVP.
# Run from repository root: Rscript tools/make_example_data.R

repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(repo_root, "R"))) {
  stop("Please run this script from the package root directory.", call. = FALSE)
}

setwd(repo_root)

source(file.path("R", "utils.R"), chdir = TRUE)
source(file.path("R", "io-seurat.R"), chdir = TRUE)
source(file.path("R", "data-build.R"), chdir = TRUE)

marker_path <- file.path("inst", "extdata", "marker_list.rds")
if (!file.exists(marker_path)) {
  dir.create(dirname(marker_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(aegis_default_markers(), marker_path)
  message("Wrote marker list: ", marker_path)
} else {
  message("Marker list exists: ", marker_path)
}

build_aegis_example_data(
  raw_dir = ".",
  output_file = file.path("data", "aegis_example.rda"),
  object_name = "aegis_example",
  max_spots = 1200L,
  seed = 42
)

message("Wrote example object: data/aegis_example.rda")
