# Plot spot-level consensus confidence map

Plot spot-level consensus confidence map

## Usage

``` r
plot_consensus_confidence(x, sample = NULL, palette = "nature", base_size = 12)
```

## Arguments

- x:

  An `aegis` object.

- sample:

  Optional sample ID for multi-sample objects.

- palette:

  Palette family: `nature`, `viridis`, `scico`, or `brewer`.

- base_size:

  Base font size.

## Value

A plot object (Seurat spatial plot when available, otherwise ggplot).
