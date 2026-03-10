# Read RCTD deconvolution outputs

Reads exported RCTD result files and standardizes them to a
spot-by-celltype numeric matrix usable by
[`as_aegis()`](https://jameswu7.github.io/AEGIS/reference/as_aegis.md).

## Usage

``` r
read_rctd(
  path,
  type = c("auto", "weights", "proportions", "results_df", "rds"),
  spot_col = NULL,
  normalize = TRUE,
  strict = TRUE
)
```

## Arguments

- path:

  File path to an exported RCTD result (csv/tsv/txt/rds).

- type:

  Input type: `auto`, `weights`, `proportions`, `results_df`, `rds`.

- spot_col:

  Optional explicit spot/barcode column name.

- normalize:

  If `TRUE`, rows are normalized to sum 1 when possible.

- strict:

  If `TRUE`, ambiguous parsing/transposition raises clear errors.

## Value

A numeric spot-by-celltype matrix.
