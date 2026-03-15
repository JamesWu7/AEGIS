# Plot method comparison and consensus outputs

Plot method comparison and consensus outputs

## Usage

``` r
plot_compare(
  x,
  type = c("heatmap", "spot_agreement", "consensus_map", "disagreement_map",
    "confidence_map", "ranking"),
  sample = NULL,
  palette = "nature",
  base_size = 12
)
```

## Arguments

- x:

  An `aegis` object.

- type:

  One of `heatmap`, `spot_agreement`, `consensus_map`,
  `disagreement_map`, `confidence_map`, or `ranking`.

- sample:

  Optional sample ID for multi-sample objects.

- palette:

  Palette family: `nature`, `viridis`, `scico`, or `brewer`.

- base_size:

  Base font size for the plot theme.

## Value

A ggplot object.
