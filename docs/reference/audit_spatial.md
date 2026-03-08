# Run spatial smoothness audit

Uses nearest-neighbor neighborhoods on tissue coordinates to estimate
local inconsistency of predicted cell-type proportions.

## Usage

``` r
audit_spatial(x, k = 6L)
```

## Arguments

- x:

  An `aegis` object.

- k:

  Number of nearest neighbors.

## Value

The modified `aegis` object with `x$audit$spatial` populated.
