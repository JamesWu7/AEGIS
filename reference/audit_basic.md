# Run basic quality audits on deconvolution matrices

Computes per-spot diagnostics and method-level sparsity summaries.

## Usage

``` r
audit_basic(x, threshold = 0.05)
```

## Arguments

- x:

  An `aegis` object.

- threshold:

  Numeric threshold for counting detected cell types.

## Value

The modified `aegis` object with `x$audit$basic` populated.
