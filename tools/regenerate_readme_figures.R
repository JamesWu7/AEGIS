#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
if (length(file_arg) > 0L) {
  script_path <- sub("^--file=", "", file_arg[[1L]])
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(".", mustWork = TRUE)
}
setwd(repo_root)

required_pkgs <- c("devtools", "ggplot2", "Seurat")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop(sprintf("Missing required packages for figure regeneration: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

devtools::load_all(".", quiet = TRUE)

data("aegis_example", package = "AEGIS")
seu <- aegis_example
markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))
deconv <- simulate_deconv_results(seu, seed = 2026)
obj <- run_aegis(seu, deconv = deconv, markers = markers)

out_dir <- file.path("inst", "assets", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Slice context figure.
p_slice <- Seurat::SpatialFeaturePlot(seu, features = "nCount_Spatial")
ggplot2::ggsave(
  filename = file.path(out_dir, "readme-slice.png"),
  plot = p_slice,
  width = 7.5,
  height = 6,
  dpi = 180
)

# Dominance figure from package plotting API.
p_dom <- plot_audit(obj, type = "dominance", method = "RCTD", palette = "nature")
ggplot2::ggsave(
  filename = file.path(out_dir, "readme-dominance.png"),
  plot = p_dom,
  width = 7.5,
  height = 5.5,
  dpi = 180
)

# Comparison heatmap.
p_heat <- plot_compare(obj, type = "heatmap", palette = "nature")
ggplot2::ggsave(
  filename = file.path(out_dir, "readme-heatmap.png"),
  plot = p_heat,
  width = 8,
  height = 5.5,
  dpi = 180
)

# Consensus map.
p_cons <- plot_compare(obj, type = "consensus_map", palette = "brewer")
ggplot2::ggsave(
  filename = file.path(out_dir, "readme-consensus.png"),
  plot = p_cons,
  width = 7.5,
  height = 5.5,
  dpi = 180
)

message("README figures regenerated under inst/assets/figures/")
