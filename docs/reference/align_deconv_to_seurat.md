# Align a deconvolution matrix to Seurat spot order

Align a deconvolution matrix to Seurat spot order

## Usage

``` r
align_deconv_to_seurat(deconv, seu, method_name = "deconv")
```

## Arguments

- deconv:

  A deconvolution matrix/data.frame with spot row names.

- seu:

  A Seurat object.

- method_name:

  Method label used in error/warning messages.

## Value

A numeric matrix aligned to `colnames(seu)`.
