#!/usr/bin/env Rscript

repo_root <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(repo_root, "R"))) {
  stop("Run this script from package root.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(AEGIS)
  library(ggplot2)
  library(Seurat)
})

data("aegis_example", package = "AEGIS")
markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))
all_methods <- get_supported_methods()$method_name

deconv <- simulate_deconv_results(
  aegis_example,
  methods = all_methods,
  cell_types = c("B_cell", "T_cell", "Myeloid"),
  seed = 123
)
obj <- as_aegis(aegis_example, deconv, markers = markers)
obj <- audit_basic(obj)
obj <- audit_marker(obj)
obj <- audit_spatial(obj)
obj <- compare_methods(obj)
obj <- score_methods(obj)
obj <- rank_methods(obj, method = "mean_rank")
obj <- compute_consensus(obj, strategy = "weighted", top_n = min(4, length(all_methods)))

out_dir <- file.path("inst", "assets", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Slice-level overview plot for the spatial transcriptomics section.
p_slice <- Seurat::SpatialFeaturePlot(aegis_example, features = "nCount_Spatial") +
  ggplot2::ggtitle("Human Lymph Node slice (nCount_Spatial)")

ggsave(file.path(out_dir, "readme-slice.png"), p_slice, width = 7.2, height = 6.2, dpi = 160)
ggsave(file.path(out_dir, "readme-dominance.png"), plot_audit(obj, "dominance"), width = 9, height = 5, dpi = 160)
ggsave(file.path(out_dir, "readme-heatmap.png"), plot_compare(obj, "heatmap"), width = 15, height = 6.5, dpi = 170)
ggsave(file.path(out_dir, "readme-ranking.png"), plot_compare(obj, "ranking"), width = 9, height = 6.2, dpi = 170)
ggsave(file.path(out_dir, "readme-spot-agreement.png"), plot_compare(obj, "spot_agreement"), width = 7.2, height = 5.4, dpi = 160)
ggsave(file.path(out_dir, "readme-confidence.png"), plot_compare(obj, "confidence_map"), width = 7.2, height = 5.4, dpi = 160)

message("README figures generated at: ", normalizePath(out_dir))
