# Construct a multi-sample AEGIS object

Construct a multi-sample AEGIS object

## Usage

``` r
as_aegis_multi(seu_list, deconv, markers = NULL, meta = NULL, strict = TRUE)
```

## Arguments

- seu_list:

  Named list of Seurat spatial objects.

- deconv:

  Nested named list of deconvolution matrices by sample.

- markers:

  Optional marker list.

- meta:

  Optional metadata list.

- strict:

  If `TRUE`, enforce strict validation.

## Value

An object of class `aegis_multi`.
