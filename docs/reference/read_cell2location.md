# Read cell2location deconvolution outputs

Read cell2location deconvolution outputs

## Usage

``` r
read_cell2location(path, spot_col = NULL, normalize = TRUE, strict = TRUE)
```

## Arguments

- path:

  File path to an exported cell2location result (csv/tsv/txt/rds).

- spot_col:

  Optional explicit spot/barcode column name.

- normalize:

  If `TRUE`, abundances are row-normalized to proportions.

- strict:

  If `TRUE`, ambiguous parsing/transposition raises clear errors.

## Value

A numeric spot-by-celltype matrix.
