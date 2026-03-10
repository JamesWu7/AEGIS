# AEGIS

![AEGIS logo](inst/assets/AEGIS_Logo.jpg)

**AEGIS**: **A**udit and **E**valuate deconvolution outputs in
**G**rid-based Spatial transcriptomics.

[![R-CMD-check](https://github.com/JamesWu7/AEGIS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JamesWu7/AEGIS/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/JamesWu7/AEGIS/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/JamesWu7/AEGIS/actions/workflows/pkgdown.yaml)

AEGIS is an R package for basic auditing of spatial deconvolution
outputs on Seurat spatial objects, with a minimal and reproducible Human
Lymph Node workflow.

## Installation

``` r
install.packages("devtools")
devtools::install_github("JamesWu7/AEGIS")
```

## Workflow at a Glance

AEGIS now supports two primary workflows:

1.  Simulated method outputs for development and demos
    ([`simulate_deconv_results()`](https://jameswu7.github.io/AEGIS/reference/simulate_deconv_results.md)).
2.  Real exported outputs from external methods
    ([`read_rctd()`](https://jameswu7.github.io/AEGIS/reference/read_rctd.md),
    [`read_spotlight()`](https://jameswu7.github.io/AEGIS/reference/read_spotlight.md),
    [`read_cell2location()`](https://jameswu7.github.io/AEGIS/reference/read_cell2location.md)).

For day-to-day use, the recommended minimal API is:

1.  [`load_10x_lymphnode()`](https://jameswu7.github.io/AEGIS/reference/load_10x_lymphnode.md)
    or
    [`load_10x_spatial_set()`](https://jameswu7.github.io/AEGIS/reference/load_10x_spatial_set.md)
2.  [`run_aegis()`](https://jameswu7.github.io/AEGIS/reference/run_aegis.md)
3.  [`plot_audit()`](https://jameswu7.github.io/AEGIS/reference/plot_audit.md)
    /
    [`plot_compare()`](https://jameswu7.github.io/AEGIS/reference/plot_compare.md)
    /
    [`render_report()`](https://jameswu7.github.io/AEGIS/reference/render_report.md)

## Quick Start (Simulated)

``` r
library(AEGIS)

seu <- load_10x_lymphnode()
deconv <- simulate_deconv_results(seu)
markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))
obj <- run_aegis(seu, deconv = deconv, markers = markers)
```

## Import Real Deconvolution Results (P5)

AEGIS imports exported result tables from external methods. It does
**not** install or run RCTD/SPOTlight/cell2location backends.

``` r
seu <- load_10x_lymphnode()

rctd <- read_rctd("path/to/rctd_output.csv")
spotlight <- read_spotlight("path/to/spotlight_output.tsv")
cell2location <- read_cell2location("path/to/cell2location_output.csv")

obj <- as_aegis(
  seu,
  deconv = list(
    RCTD = rctd,
    SPOTlight = spotlight,
    cell2location = cell2location
  )
)

obj <- audit_basic(obj)
obj <- compare_methods(obj)
obj <- compute_consensus(obj)
```

For cell2location, export posterior abundance/proportion tables to
csv/tsv/txt first, then import with
[`read_cell2location()`](https://jameswu7.github.io/AEGIS/reference/read_cell2location.md).

## Multi-sample Workflow (P6)

``` r
seu_list <- load_10x_spatial_set(
  paths = c("sample1_dir", "sample2_dir"),
  sample_ids = c("sample1", "sample2")
)

deconv_nested <- list(
  sample1 = list(RCTD = rctd1, SPOTlight = spotlight1),
  sample2 = list(RCTD = rctd2, SPOTlight = spotlight2)
)

obj_multi <- run_aegis(seu_list, deconv = deconv_nested, markers = markers)

summary_tbl <- summarize_by_sample(obj_multi)
render_report_batch(obj_multi, output_dir = "reports")
```

## Complete Tutorials

If GitHub Pages is temporarily unavailable, use the preview fallback
links or the source `.Rmd` links below.

- [Overview tutorial (object model +
  workflows)](https://jameswu7.github.io/AEGIS/articles/AEGIS-overview.html)
  ([preview
  fallback](https://htmlpreview.github.io/?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-overview.html),
  [source](https://jameswu7.github.io/AEGIS/vignettes/AEGIS-overview.Rmd))
- [Human lymph node demo (end-to-end
  demo)](https://jameswu7.github.io/AEGIS/articles/AEGIS-demo-human-lymph-node.html)
  ([preview
  fallback](https://htmlpreview.github.io/?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-demo-human-lymph-node.html),
  [source](https://jameswu7.github.io/AEGIS/vignettes/AEGIS-demo-human-lymph-node.Rmd))
- [Complete tutorial (simulated + real import +
  multi-sample)](https://jameswu7.github.io/AEGIS/articles/AEGIS-complete-tutorial.html)
  ([preview
  fallback](https://htmlpreview.github.io/?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-complete-tutorial.html),
  [source](https://jameswu7.github.io/AEGIS/vignettes/AEGIS-complete-tutorial.Rmd))

## Key Functions

- [`load_10x_lymphnode()`](https://jameswu7.github.io/AEGIS/reference/load_10x_lymphnode.md):
  load the Human Lymph Node 10x spatial dataset into a Seurat object.
- [`simulate_deconv_results()`](https://jameswu7.github.io/AEGIS/reference/simulate_deconv_results.md):
  generate realistic mock method outputs (spot-by-celltype proportions).
- [`read_rctd()`](https://jameswu7.github.io/AEGIS/reference/read_rctd.md):
  import exported RCTD result tables/RDS and standardize to
  spot-by-celltype.
- [`read_spotlight()`](https://jameswu7.github.io/AEGIS/reference/read_spotlight.md):
  import exported SPOTlight result tables/RDS and standardize to
  spot-by-celltype.
- [`read_cell2location()`](https://jameswu7.github.io/AEGIS/reference/read_cell2location.md):
  import exported cell2location tables/RDS (abundance or proportion) and
  standardize.
- [`as_aegis()`](https://jameswu7.github.io/AEGIS/reference/as_aegis.md):
  validate inputs and create the internal `aegis` S3 object.
- [`audit_basic()`](https://jameswu7.github.io/AEGIS/reference/audit_basic.md):
  compute per-spot and per-method basic quality metrics.
- [`audit_marker()`](https://jameswu7.github.io/AEGIS/reference/audit_marker.md):
  quantify marker-expression support and method concordance.
- [`audit_spatial()`](https://jameswu7.github.io/AEGIS/reference/audit_spatial.md):
  compute neighborhood-based local inconsistency metrics.
- [`compare_methods()`](https://jameswu7.github.io/AEGIS/reference/compare_methods.md):
  summarize cross-method agreement by cell type and spot.
- [`compute_consensus()`](https://jameswu7.github.io/AEGIS/reference/compute_consensus.md):
  aggregate shared cell types and derive confidence/stability.
- [`run_aegis()`](https://jameswu7.github.io/AEGIS/reference/run_aegis.md):
  one-call pipeline for single-sample or multi-sample workflows.

## Example Figures

### Spatial transcriptomics slice (Human Lymph Node)

![Human Lymph Node slice](inst/assets/figures/readme-slice.png)

Human Lymph Node slice

### Basic audit (dominance)

![Dominance spatial map](inst/assets/figures/readme-dominance.png)

Dominance spatial map

## Citation

``` r
citation("AEGIS")
```

BibTeX:

``` bibtex
@Manual{Wu2026AEGIS,
  title = {AEGIS: Audit and Evaluate Deconvolution Outputs in Grid-Based Spatial Transcriptomics},
  author = {Xinjie Wu},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/JamesWu7/AEGIS}
}
```
