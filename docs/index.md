# AEGIS

![AEGIS logo](inst/assets/AEGIS_Logo.jpg)

**AEGIS**: **A**udit and **E**valuate deconvolution outputs in
**G**rid-based Spatial transcriptomics.

[![R-CMD-check](https://github.com/JamesWu7/AEGIS/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/JamesWu7/AEGIS/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/JamesWu7/AEGIS/actions/workflows/pkgdown.yaml/badge.svg?branch=main)](https://github.com/JamesWu7/AEGIS/actions/workflows/pkgdown.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

AEGIS is an R package for basic auditing of spatial deconvolution
outputs on Seurat spatial objects, with a minimal and reproducible Human
Lymph Node workflow.

## Installation

``` r
install.packages("devtools")
devtools::install_github("JamesWu7/AEGIS")
```

## Workflow at a Glance

AEGIS now supports three primary workflows:

1.  Simulated method outputs for development and demos
    ([`simulate_deconv_results()`](https://jameswu7.github.io/AEGIS/reference/simulate_deconv_results.md)).
2.  Real exported outputs from external methods through adapter readers
    (e.g. [`read_rctd()`](https://jameswu7.github.io/AEGIS/reference/read_rctd.md),
    [`read_spotlight()`](https://jameswu7.github.io/AEGIS/reference/read_spotlight.md),
    [`read_cell2location()`](https://jameswu7.github.io/AEGIS/reference/read_cell2location.md),
    [`read_card()`](https://jameswu7.github.io/AEGIS/reference/read_card.md),
    [`read_destvi()`](https://jameswu7.github.io/AEGIS/reference/read_destvi.md),
    [`read_stdeconvolve()`](https://jameswu7.github.io/AEGIS/reference/read_stdeconvolve.md)),
    plus
    [`read_deconv_table()`](https://jameswu7.github.io/AEGIS/reference/read_deconv_table.md)
    for generic spot-by-celltype tables.
3.  One-shot deconvolution orchestration for directly runnable methods
    through
    [`run_deconvolution()`](https://jameswu7.github.io/AEGIS/reference/run_deconvolution.md)
    and
    [`run_aegis_full()`](https://jameswu7.github.io/AEGIS/reference/run_deconvolution.md),
    with capability metadata from
    [`get_supported_methods()`](https://jameswu7.github.io/AEGIS/reference/get_supported_methods.md).

For day-to-day use, the recommended minimal API is:

1.  [`load_10x_lymphnode()`](https://jameswu7.github.io/AEGIS/reference/load_10x_lymphnode.md)
    or
    [`load_10x_spatial_set()`](https://jameswu7.github.io/AEGIS/reference/load_10x_spatial_set.md)
2.  [`run_aegis()`](https://jameswu7.github.io/AEGIS/reference/run_aegis.md)
3.  [`score_methods()`](https://jameswu7.github.io/AEGIS/reference/score_methods.md)
    -\>
    [`rank_methods()`](https://jameswu7.github.io/AEGIS/reference/rank_methods.md)
    -\> `compute_consensus(strategy = "weighted")`
4.  [`plot_method_ranking()`](https://jameswu7.github.io/AEGIS/reference/plot_method_ranking.md)
    /
    [`plot_disagreement_map()`](https://jameswu7.github.io/AEGIS/reference/plot_disagreement_map.md)
    /
    [`plot_consensus_confidence()`](https://jameswu7.github.io/AEGIS/reference/plot_consensus_confidence.md)
    /
    [`render_report()`](https://jameswu7.github.io/AEGIS/reference/render_report.md)

## Quick Start (Simulated)

``` r
library(AEGIS)

seu <- load_10x_lymphnode()
deconv <- simulate_deconv_results(seu)
markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))
obj <- run_aegis(seu, deconv = deconv, markers = markers)
obj <- score_methods(obj)
obj <- rank_methods(obj, method = "rra")
obj <- compute_consensus(obj, strategy = "weighted")

p_rank <- plot_method_ranking(obj)
p_dis <- plot_disagreement_map(obj)
p_conf <- plot_consensus_confidence(obj)
```

## One-shot Deconvolution (P9)

``` r
seu <- load_10x_lymphnode()

res <- run_deconvolution(
  seu = seu,
  reference = ref,
  methods = c("SPOTlight", "RCTD", "CARD"),
  strict = FALSE
)

obj <- run_aegis(res$seu, deconv = res$deconv, markers = markers)
```

Use
[`get_supported_methods()`](https://jameswu7.github.io/AEGIS/reference/get_supported_methods.md)
to inspect support mode (`run_and_import_r`, `run_and_import_python`,
`import_only`) before execution.

## Import Real Deconvolution Results (P8)

AEGIS imports exported result tables from external methods. It does
**not** install or run external backends. For Python/deep-learning
methods, export spot-by-celltype tables first, then import into AEGIS.

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

### Method Support Matrix

| Method | Support mode | Run in R | Run via Python | Import exported results | Notes |
|---|---|---|---|---|---|
| RCTD | `run_and_import_r` | Yes (`run_rctd`) | No | Yes (`read_rctd`) | Direct run requires `spacexr` |
| SPOTlight | `run_and_import_r` | Yes (`run_spotlight`) | No | Yes (`read_spotlight`) | Direct run requires `SPOTlight` |
| CARD | `run_and_import_r` | Yes (`run_card`) | No | Yes (`read_card`) | Direct run requires `CARD` |
| cell2location | `run_and_import_python` | No | Optional (`run_cell2location`) | Yes (`read_cell2location`) | Python/reticulate environment required |
| stereoscope | `run_and_import_python` | No | Optional (`run_stereoscope`) | Yes (`read_stereoscope`) | Python/reticulate environment required |
| DestVI | `run_and_import_python` | No | Optional (`run_destvi`) | Yes (`read_destvi`) | Python/reticulate environment required |
| Tangram | `run_and_import_python` | No | Optional (`run_tangram`) | Yes (`read_tangram`) | Mapping-style composition input |
| SpatialDWLS | `import_only` | No | No | Yes (`read_spatialdwls`) | Table adapter |
| STdeconvolve | `import_only` | No | No | Yes (`read_stdeconvolve`) | Latent labels allowed |
| DSTG | `import_only` | No | No | Yes (`read_dstg`) | Table adapter |
| STRIDE | `import_only` | No | No | Yes (`read_stride`) | Topic-only strict checks supported |

### Additional Import Examples

``` r
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

## Tutorials

If GitHub Pages is temporarily unavailable, use the preview fallback
links or the source `.Rmd` links below.

- [Quick Start
  tutorial](https://htmlpreview.github.io/?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-overview.html)
  (includes `plot_compare` visualizations and RRA/mean-rank selection)
  ([pkgdown
  page](https://jameswu7.github.io/AEGIS/articles/AEGIS-overview.html),
  [source](https://jameswu7.github.io/AEGIS/vignettes/AEGIS-overview.Rmd))
- [One-step deconvolution
  tutorial](https://htmlpreview.github.io/?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-one-step-deconvolution.html)
  (focuses on `get_supported_methods()`, `run_deconvolution()`,
  `run_aegis_full()`, and practical run-vs-import decisions) ([pkgdown
  page](https://jameswu7.github.io/AEGIS/articles/AEGIS-one-step-deconvolution.html),
  [source](https://jameswu7.github.io/AEGIS/vignettes/AEGIS-one-step-deconvolution.Rmd))
- [Deconvolution from Scratch + Real Data tutorial (Human Lymph
  Node)](https://htmlpreview.github.io/?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-complete-tutorial.html)
  (includes `get_supported_methods()`, `run_deconvolution()`,
  `run_aegis_full()`, all supported import adapters, method ranking,
  weighted consensus, and `plot_compare`-based visualization) ([pkgdown
  page](https://jameswu7.github.io/AEGIS/articles/AEGIS-complete-tutorial.html),
  [source](https://jameswu7.github.io/AEGIS/vignettes/AEGIS-complete-tutorial.Rmd))

Local preview HTML files are also kept in `vignettes/` for direct
viewing: `vignettes/AEGIS-overview.html`,
`vignettes/AEGIS-one-step-deconvolution.html`, and
`vignettes/AEGIS-complete-tutorial.html`.

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
- [`read_card()`](https://jameswu7.github.io/AEGIS/reference/read_card.md),
  [`read_spatialdwls()`](https://jameswu7.github.io/AEGIS/reference/read_spatialdwls.md),
  [`read_stereoscope()`](https://jameswu7.github.io/AEGIS/reference/read_stereoscope.md),
  [`read_destvi()`](https://jameswu7.github.io/AEGIS/reference/read_destvi.md),
  [`read_tangram()`](https://jameswu7.github.io/AEGIS/reference/read_tangram.md),
  [`read_stdeconvolve()`](https://jameswu7.github.io/AEGIS/reference/read_stdeconvolve.md),
  [`read_dstg()`](https://jameswu7.github.io/AEGIS/reference/read_dstg.md),
  [`read_stride()`](https://jameswu7.github.io/AEGIS/reference/read_stride.md):
  method-specific import adapters.
- [`read_deconv_table()`](https://jameswu7.github.io/AEGIS/reference/read_deconv_table.md):
  generic importer for spot-by-celltype exported tables.
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
- [`score_methods()`](https://jameswu7.github.io/AEGIS/reference/score_methods.md):
  score methods from marker/spatial/agreement/stability evidence.
- [`rank_methods()`](https://jameswu7.github.io/AEGIS/reference/rank_methods.md):
  aggregate evidence into robust method rankings (`rra` or `mean_rank`).
- [`compute_consensus()`](https://jameswu7.github.io/AEGIS/reference/compute_consensus.md):
  integrate methods with `mean` / `weighted` / `trimmed_mean` strategies
  and return disagreement/confidence.
- [`plot_method_ranking()`](https://jameswu7.github.io/AEGIS/reference/plot_method_ranking.md):
  ggplot ranking summary from strongest to weakest methods.
- [`plot_disagreement_map()`](https://jameswu7.github.io/AEGIS/reference/plot_disagreement_map.md):
  tissue-context map of spot-level cross-method disagreement.
- [`plot_consensus_confidence()`](https://jameswu7.github.io/AEGIS/reference/plot_consensus_confidence.md):
  tissue-context map of spot-level consensus confidence.
- [`run_aegis()`](https://jameswu7.github.io/AEGIS/reference/run_aegis.md):
  one-call pipeline for single-sample or multi-sample workflows.

## Example Figures

### Spatial transcriptomics slice (Human Lymph Node)

![Human Lymph Node slice](inst/assets/figures/readme-slice.png)

Human Lymph Node slice

### Basic audit (dominance)

![Dominance spatial map](inst/assets/figures/readme-dominance.png)

Dominance spatial map

### Cross-method agreement heatmap

![Method agreement heatmap](inst/assets/figures/readme-heatmap.png)

Method agreement heatmap

### Method ranking

![Method ranking](inst/assets/figures/readme-ranking.png)

Method ranking

### Spot-level agreement

![Spot-level agreement](inst/assets/figures/readme-spot-agreement.png)

Spot-level agreement

### Consensus confidence

![Consensus confidence](inst/assets/figures/readme-confidence.png)

Consensus confidence

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
