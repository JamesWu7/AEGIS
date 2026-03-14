<p align="center">
  <img src="inst/assets/AEGIS_Logo.jpg" alt="AEGIS logo" width="360"/>
</p>

# AEGIS

**AEGIS**: **A**udit and **E**valuate deconvolution outputs in **G**rid-based Spatial transcriptomics.

[![R-CMD-check](https://github.com/JamesWu7/AEGIS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JamesWu7/AEGIS/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/JamesWu7/AEGIS/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/JamesWu7/AEGIS/actions/workflows/pkgdown.yaml)

AEGIS is an R package for basic auditing of spatial deconvolution outputs on Seurat spatial objects, with a minimal and reproducible Human Lymph Node workflow.

## Installation

```r
install.packages("devtools")
devtools::install_github("JamesWu7/AEGIS")
```

## Workflow at a Glance

AEGIS now supports two primary workflows:

1. Simulated method outputs for development and demos (`simulate_deconv_results()`).
2. Real exported outputs from external methods through adapter readers (e.g. `read_rctd()`, `read_spotlight()`, `read_cell2location()`, `read_card()`, `read_destvi()`, `read_stdeconvolve()`), plus `read_deconv_table()` for generic spot-by-celltype tables.

For day-to-day use, the recommended minimal API is:

1. `load_10x_lymphnode()` or `load_10x_spatial_set()`
2. `run_aegis()`
3. `plot_audit()` / `plot_compare()` / `render_report()`

## Quick Start (Simulated)

```r
library(AEGIS)

seu <- load_10x_lymphnode()
deconv <- simulate_deconv_results(seu)
markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))
obj <- run_aegis(seu, deconv = deconv, markers = markers)
```

## Import Real Deconvolution Results (P8)

AEGIS imports exported result tables from external methods. It does **not** install or run external backends.
For Python/deep-learning methods, export spot-by-celltype tables first, then import into AEGIS.

```r
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

For cell2location, export posterior abundance/proportion tables to csv/tsv/txt first, then import with `read_cell2location()`.

### Method Support Matrix

| Method | Adapter | Expected input | `normalize` | Notes |
|---|---|---|---|---|
| RCTD | `read_rctd()` | csv/tsv/txt/rds | Yes | Supports common exported table/RDS forms |
| SPOTlight | `read_spotlight()` | csv/tsv/txt/rds | Yes | Drops obvious metadata columns |
| cell2location | `read_cell2location()` | csv/tsv/txt/rds | Yes | Abundance tables supported |
| CARD | `read_card()` | csv/tsv/txt/rds | Yes | Table adapter |
| SpatialDWLS | `read_spatialdwls()` | csv/tsv/txt/rds | Yes | Table adapter |
| stereoscope | `read_stereoscope()` | csv/tsv/txt/rds | Yes | Export table only |
| DestVI | `read_destvi()` | csv/tsv/txt/rds | Yes | Export table only |
| Tangram | `read_tangram()` | csv/tsv/txt/rds | Yes | Treated as mapping-derived composition input |
| STdeconvolve | `read_stdeconvolve()` | csv/tsv/txt/rds | Yes | Latent labels (e.g., `topic1`) allowed |
| DSTG | `read_dstg()` | csv/tsv/txt/rds | Yes | Export table only |
| STRIDE | `read_stride()` | csv/tsv/txt/rds | Yes | Topic-only files fail in strict mode |

### Additional Import Examples

```r
card <- read_card("path/to/card.csv")
destvi <- read_destvi("path/to/destvi.csv")
stdec <- read_stdeconvolve("path/to/stdeconvolve.csv")

obj <- as_aegis(
  seu,
  deconv = list(
    CARD = card,
    DestVI = destvi,
    STdeconvolve = stdec
  )
)
```

## Multi-sample Workflow (P6)

```r
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

If GitHub Pages is temporarily unavailable, use the preview fallback links or the source `.Rmd` links below.

- [Overview tutorial (object model + workflows)](https://htmlpreview.github.io/?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-overview.html) ([pkgdown page](https://jameswu7.github.io/AEGIS/articles/AEGIS-overview.html), [source](vignettes/AEGIS-overview.Rmd))
- [Human lymph node demo (end-to-end demo)](https://htmlpreview.github.io/?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-demo-human-lymph-node.html) ([pkgdown page](https://jameswu7.github.io/AEGIS/articles/AEGIS-demo-human-lymph-node.html), [source](vignettes/AEGIS-demo-human-lymph-node.Rmd))
- [Complete tutorial (simulated + real import + multi-sample)](https://htmlpreview.github.io/?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-complete-tutorial.html) ([pkgdown page](https://jameswu7.github.io/AEGIS/articles/AEGIS-complete-tutorial.html), [source](vignettes/AEGIS-complete-tutorial.Rmd))

## Key Functions

- `load_10x_lymphnode()`: load the Human Lymph Node 10x spatial dataset into a Seurat object.
- `simulate_deconv_results()`: generate realistic mock method outputs (spot-by-celltype proportions).
- `read_rctd()`: import exported RCTD result tables/RDS and standardize to spot-by-celltype.
- `read_spotlight()`: import exported SPOTlight result tables/RDS and standardize to spot-by-celltype.
- `read_cell2location()`: import exported cell2location tables/RDS (abundance or proportion) and standardize.
- `read_card()`, `read_spatialdwls()`, `read_stereoscope()`, `read_destvi()`, `read_tangram()`, `read_stdeconvolve()`, `read_dstg()`, `read_stride()`: method-specific import adapters.
- `read_deconv_table()`: generic importer for spot-by-celltype exported tables.
- `as_aegis()`: validate inputs and create the internal `aegis` S3 object.
- `audit_basic()`: compute per-spot and per-method basic quality metrics.
- `audit_marker()`: quantify marker-expression support and method concordance.
- `audit_spatial()`: compute neighborhood-based local inconsistency metrics.
- `compare_methods()`: summarize cross-method agreement by cell type and spot.
- `compute_consensus()`: aggregate shared cell types and derive confidence/stability.
- `run_aegis()`: one-call pipeline for single-sample or multi-sample workflows.

## Example Figures

### Spatial transcriptomics slice (Human Lymph Node)

![Human Lymph Node slice](inst/assets/figures/readme-slice.png)

### Basic audit (dominance)

![Dominance spatial map](inst/assets/figures/readme-dominance.png)

## Citation

```r
citation("AEGIS")
```

BibTeX:

```bibtex
@Manual{Wu2026AEGIS,
  title = {AEGIS: Audit and Evaluate Deconvolution Outputs in Grid-Based Spatial Transcriptomics},
  author = {Xinjie Wu},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/JamesWu7/AEGIS}
}
```
