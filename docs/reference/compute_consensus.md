# Compute consensus deconvolution profile

Aggregates deconvolution outputs across methods (mean aggregation for
MVP) and computes disagreement/stability summaries.

## Usage

``` r
compute_consensus(x)
```

## Arguments

- x:

  An `aegis` object.

## Value

Modified `aegis` object with `x$consensus$result` populated.
