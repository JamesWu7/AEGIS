# Simulate deconvolution results for a Seurat spatial object

Generates realistic mock spot-by-celltype proportion matrices for
multiple methods. This is intended for MVP development and demos before
plugging in real deconvolution backends.

## Usage

``` r
simulate_deconv_results(
  seu,
  cell_types = c("B_cell", "T_cell", "NK_cell", "Plasma", "Myeloid", "Endothelial",
    "Stromal"),
  methods = c("RCTD", "SPOTlight", "cell2location"),
  seed = 123,
  noise_scale = 0.05
)
```

## Arguments

- seu:

  A Seurat object.

- cell_types:

  Character vector of cell types.

- methods:

  Character vector of method names.

- seed:

  Random seed for reproducibility.

- noise_scale:

  Numeric scale for method-specific perturbation.

## Value

A named list of numeric matrices (spot-by-celltype proportions).
