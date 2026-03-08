<p align="center">
  <img src="inst/assets/AEGIS_Logo.jpg" alt="AEGIS logo" width="360"/>
</p>

# AEGIS

AEGIS is an R package for basic auditing of spatial deconvolution outputs on Seurat spatial objects, with a minimal and reproducible Human Lymph Node workflow.

## Installation

```r
install.packages("devtools")
devtools::install_github("JamesWu7/AEGIS")
```

## Quick Start

```r
library(AEGIS)

seu <- load_10x_lymphnode()
deconv <- simulate_deconv_results(seu)
obj <- as_aegis(seu, deconv)
obj <- audit_basic(obj)
```

## Complete Tutorials

- [Overview tutorial](https://github.com/JamesWu7/AEGIS/blob/main/vignettes/AEGIS-overview.Rmd)
- [Human lymph node demo](https://github.com/JamesWu7/AEGIS/blob/main/vignettes/AEGIS-demo-human-lymph-node.Rmd)
- [Complete tutorial](https://github.com/JamesWu7/AEGIS/blob/main/vignettes/AEGIS-complete-tutorial.Rmd)

## Key Functions

- `load_10x_lymphnode()`: load the Human Lymph Node 10x spatial dataset into a Seurat object.
- `simulate_deconv_results()`: generate realistic mock method outputs (spot-by-celltype proportions).
- `as_aegis()`: validate inputs and create the internal `aegis` S3 object.
- `audit_basic()`: compute per-spot and per-method basic quality metrics.

## Example Figures

### Spatial transcriptomics slice (Human Lymph Node)

![Human Lymph Node slice](inst/assets/figures/readme-slice.png)

### Basic audit (dominance)

![Dominance spatial map](inst/assets/figures/readme-dominance.png)
