# Load a set of 10x spatial datasets

Loads multiple spatial datasets into a named list of Seurat objects and
attaches sample-aware metadata.

## Usage

``` r
load_10x_spatial_set(
  paths,
  sample_ids = NULL,
  section_ids = NULL,
  strict = TRUE
)
```

## Arguments

- paths:

  Character vector of sample directories.

- sample_ids:

  Optional sample identifiers.

- section_ids:

  Optional section identifiers.

- strict:

  If `TRUE`, ambiguous ID inference errors clearly.

## Value

Named list of Seurat objects.
