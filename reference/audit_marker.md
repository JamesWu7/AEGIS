# Run marker-based concordance audit

Computes marker expression scores per cell type and compares them with
deconvolution abundance estimates.

## Usage

``` r
audit_marker(x, markers = x$markers, assay = "Spatial")
```

## Arguments

- x:

  An `aegis` object.

- markers:

  Named list of marker genes by cell type. Defaults to `x$markers`.

- assay:

  Assay name to use from Seurat.

## Value

The modified `aegis` object with `x$audit$marker` populated.
