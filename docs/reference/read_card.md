# Read CARD deconvolution outputs

Read CARD deconvolution outputs

## Usage

``` r
read_card(path, spot_col = NULL, normalize = TRUE, strict = TRUE)
```

## Arguments

- path:

  File path to exported CARD results (csv/tsv/txt/rds).

- spot_col:

  Optional explicit spot/barcode column name.

- normalize:

  If `TRUE`, rows are normalized to sum 1 when possible.

- strict:

  If `TRUE`, ambiguous parsing/transposition raises clear errors.

## Value

A numeric spot-by-celltype matrix.
