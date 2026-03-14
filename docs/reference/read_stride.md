# Read STRIDE deconvolution outputs

Read STRIDE deconvolution outputs

## Usage

``` r
read_stride(path, spot_col = NULL, normalize = TRUE, strict = TRUE)
```

## Arguments

- path:

  File path to exported STRIDE results (csv/tsv/txt/rds).

- spot_col:

  Optional explicit spot/barcode column name.

- normalize:

  If `TRUE`, rows are normalized to sum 1 when possible.

- strict:

  If `TRUE`, ambiguous topic-only files raise clear errors.

## Value

A numeric spot-by-celltype matrix.
