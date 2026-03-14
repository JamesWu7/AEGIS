# Compute consensus deconvolution profile

Aggregates deconvolution outputs across methods and computes
disagreement and confidence summaries.

## Usage

``` r
compute_consensus(
  x,
  strategy = c("mean", "weighted", "trimmed_mean"),
  top_n = NULL
)
```

## Arguments

- x:

  An `aegis` object.

- strategy:

  Aggregation strategy: `mean`, `weighted`, or `trimmed_mean`.

- top_n:

  Optional number of top-ranked methods to include.

## Value

Modified `aegis` object with `x$consensus$result` populated.
