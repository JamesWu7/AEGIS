# Merge a list of spatial Seurat objects

Merge a list of spatial Seurat objects

## Usage

``` r
merge_spatial_seurat_list(
  seu_list,
  sample_ids = NULL,
  add_cell_ids = TRUE,
  preserve_images = TRUE,
  strict = TRUE
)
```

## Arguments

- seu_list:

  Named list of Seurat objects.

- sample_ids:

  Optional sample identifiers.

- add_cell_ids:

  If `TRUE`, prefixes spots with sample IDs before merge.

- preserve_images:

  If `TRUE`, checks image preservation and warns/errors.

- strict:

  If `TRUE`, enforce strict validation.

## Value

A merged Seurat object with `sample_id` metadata.
