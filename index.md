# AEGIS

![AEGIS logo](inst/assets/AEGIS_Logo.jpg)

AEGIS is an R package for basic auditing of spatial deconvolution
outputs on Seurat spatial objects, with a minimal and reproducible Human
Lymph Node workflow.

## Installation

``` r
install.packages("devtools")
devtools::install_github("JamesWu7/AEGIS")
```

## Quick Start

``` r
library(AEGIS)

seu <- load_10x_lymphnode()
deconv <- simulate_deconv_results(seu)
obj <- as_aegis(seu, deconv)
obj <- audit_basic(obj)
```

## Complete Tutorials

- [Overview
  tutorial](https://jameswu7.github.io/AEGIS/articles/AEGIS-overview.html)
- [Human lymph node
  demo](https://jameswu7.github.io/AEGIS/articles/AEGIS-demo-human-lymph-node.html)
- [Complete
  tutorial](https://jameswu7.github.io/AEGIS/articles/AEGIS-complete-tutorial.html)

## Key Functions

- [`load_10x_lymphnode()`](https://jameswu7.github.io/AEGIS/reference/load_10x_lymphnode.md):
  load the Human Lymph Node 10x spatial dataset into a Seurat object.
- [`simulate_deconv_results()`](https://jameswu7.github.io/AEGIS/reference/simulate_deconv_results.md):
  generate realistic mock method outputs (spot-by-celltype proportions).
- [`as_aegis()`](https://jameswu7.github.io/AEGIS/reference/as_aegis.md):
  validate inputs and create the internal `aegis` S3 object.
- [`audit_basic()`](https://jameswu7.github.io/AEGIS/reference/audit_basic.md):
  compute per-spot and per-method basic quality metrics.

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
  title = {AEGIS: Audit and Evaluate deconvolution outputs in Grid-based Spatial transcriptomics},
  author = {Xinjie Wu},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/JamesWu7/AEGIS}
}
```
