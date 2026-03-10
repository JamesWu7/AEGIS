# Run the standard AEGIS workflow with one function

Convenience entry point that constructs an AEGIS object (if needed) and
runs the core audit/comparison pipeline for single-sample and
multi-sample workflows.

## Usage

``` r
run_aegis(
  x,
  deconv = NULL,
  markers = NULL,
  do_marker = TRUE,
  do_spatial = TRUE,
  do_compare = TRUE,
  do_consensus = TRUE,
  strict = TRUE
)
```

## Arguments

- x:

  A Seurat object, named list of Seurat objects, `aegis`, or
  `aegis_multi`.

- deconv:

  Optional deconvolution input. For single-sample: named list of method
  matrices. For multi-sample: nested list by sample.

- markers:

  Optional marker list. If provided with an existing AEGIS object, it
  replaces `x$markers`.

- do_marker:

  Logical; run
  [`audit_marker()`](https://jameswu7.github.io/AEGIS/reference/audit_marker.md).

- do_spatial:

  Logical; run
  [`audit_spatial()`](https://jameswu7.github.io/AEGIS/reference/audit_spatial.md).

- do_compare:

  Logical; run
  [`compare_methods()`](https://jameswu7.github.io/AEGIS/reference/compare_methods.md).

- do_consensus:

  Logical; run
  [`compute_consensus()`](https://jameswu7.github.io/AEGIS/reference/compute_consensus.md).

- strict:

  Logical; strict validation for multi-sample constructor path.

## Value

An `aegis` or `aegis_multi` object.
