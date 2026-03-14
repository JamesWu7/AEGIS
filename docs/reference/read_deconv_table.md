# Read a generic deconvolution result table

Generic importer for spot-by-celltype tables exported by external
methods.

## Usage

``` r
read_deconv_table(
  path,
  method = NULL,
  type = c("auto", "rds", "weights", "proportions", "results_df"),
  spot_col = NULL,
  normalize = TRUE,
  strict = TRUE
)
```

## Arguments

- path:

  File path to exported results (csv/tsv/txt/rds).

- method:

  Optional method label used in messages and metadata.

- type:

  Input type for RDS-aware parsing (`auto`, `rds`, `weights`,
  `proportions`, `results_df`).

- spot_col:

  Optional explicit spot/barcode column name.

- normalize:

  If `TRUE`, rows are normalized to sum 1 when possible.

- strict:

  If `TRUE`, ambiguous parsing/transposition raises clear errors.

## Value

A numeric spot-by-celltype matrix.
