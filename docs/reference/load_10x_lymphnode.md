# Load the 10x Human Lymph Node Visium dataset as a Seurat object

This loader is tailored to the authoritative MVP assets:
`V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5`,
`V1_Human_Lymph_Node_spatial.tar.gz`, and
`V1_Human_Lymph_Node_metrics_summary.csv`.

## Usage

``` r
load_10x_lymphnode(
  data_dir = ".",
  h5_file = "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5",
  spatial_tar = "V1_Human_Lymph_Node_spatial.tar.gz",
  metrics_file = "V1_Human_Lymph_Node_metrics_summary.csv",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE
)
```

## Arguments

- data_dir:

  Directory containing the raw files.

- h5_file:

  Name of the feature-barcode H5 file.

- spatial_tar:

  Name of the spatial tar.gz archive.

- metrics_file:

  Name of the metrics CSV file.

- assay:

  Assay name passed to
  [`Seurat::Load10X_Spatial()`](https://satijalab.org/seurat/reference/Load10X_Spatial.html).

- slice:

  Slice name passed to
  [`Seurat::Load10X_Spatial()`](https://satijalab.org/seurat/reference/Load10X_Spatial.html).

- filter.matrix:

  Logical; forwarded to
  [`Seurat::Load10X_Spatial()`](https://satijalab.org/seurat/reference/Load10X_Spatial.html).

## Value

A Seurat spatial object.

## Details

It will unpack the spatial tarball automatically when `spatial/` is
missing, then call
[`Seurat::Load10X_Spatial()`](https://satijalab.org/seurat/reference/Load10X_Spatial.html).
