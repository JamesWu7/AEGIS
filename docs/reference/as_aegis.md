# Construct an AEGIS object

Creates a lightweight S3 `aegis` object that stores a Seurat object,
deconvolution matrices, marker definitions, and downstream analysis
state.

## Usage

``` r
as_aegis(seu, deconv, markers = NULL, meta = NULL)
```

## Arguments

- seu:

  A Seurat object with spatial spots in columns.

- deconv:

  Named list of deconvolution matrices/data frames. Each element must be
  spot-by-celltype with row names aligned to Seurat spot names.

- markers:

  Optional named list of marker genes by cell type.

- meta:

  Optional list of user metadata.

## Value

An object of class `aegis`.
