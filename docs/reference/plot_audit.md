# Plot audit outputs

Plot audit outputs

## Usage

``` r
plot_audit(
  x,
  type = c("sumdev", "dominance", "entropy", "marker", "smoothness"),
  method = NULL,
  sample = NULL,
  cell_type = NULL,
  palette = "nature",
  point_size = NULL,
  base_size = 12
)
```

## Arguments

- x:

  An `aegis` object.

- type:

  One of `sumdev`, `dominance`, `entropy`, `marker`, `smoothness`.

- method:

  Optional method name to subset to one deconvolution method.

- sample:

  Optional sample ID for multi-sample objects.

- cell_type:

  Optional cell type for marker/smoothness views.

- palette:

  Palette family: `nature`, `viridis`, `scico`, or `brewer`.

- point_size:

  Optional spatial point size. Auto-selected when `NULL`.

- base_size:

  Base font size for the plot theme.

## Value

A ggplot object.
