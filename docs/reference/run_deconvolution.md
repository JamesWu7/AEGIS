# Unified deconvolution orchestration and method runners

`run_deconvolution()` dispatches requested methods using the AEGIS
method registry, standardizes outputs to the shared AEGIS deconvolution
contract, and returns results ready for
[`as_aegis()`](https://jameswu7.github.io/AEGIS/reference/as_aegis.md) /
[`run_aegis()`](https://jameswu7.github.io/AEGIS/reference/run_aegis.md).

`run_aegis_full()` combines one-shot deconvolution dispatch with the
existing AEGIS downstream pipeline.

Method-specific runner wrappers are included for directly runnable
R-native methods and optional Python-backed methods.

## Usage

``` r
run_deconvolution(
  seu,
  reference,
  methods = c("SPOTlight", "RCTD", "CARD"),
  sample_id = NULL,
  normalize = TRUE,
  strict = TRUE,
  use_python = TRUE,
  runner_overrides = NULL,
  runner_args = NULL,
  ...
)

run_aegis_full(
  seu,
  reference,
  methods = c("SPOTlight", "RCTD", "CARD"),
  markers = NULL,
  strict = TRUE,
  use_python = TRUE,
  normalize = TRUE,
  ...
)

run_spotlight(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_card(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_rctd(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_cell2location(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_destvi(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_tangram(seu, reference, normalize = TRUE, strict = TRUE, ...)

run_stereoscope(seu, reference, normalize = TRUE, strict = TRUE, ...)
```

## Arguments

- seu:

  A Seurat spatial object.

- reference:

  Reference input used by runnable methods.

- methods:

  Character vector of method names.

- sample_id:

  Optional sample identifier recorded in deconvolution output metadata.

- normalize:

  Logical; row-normalize outputs to proportions when possible.

- strict:

  Logical; strict error behavior for unsupported/missing dependencies or
  failed execution.

- use_python:

  Logical; allow Python-backed methods via `reticulate`.

- runner_overrides:

  Optional named list of runner functions keyed by method, mainly for
  advanced usage/testing.

- runner_args:

  Optional named list of per-method argument lists passed to the
  selected runner.

- markers:

  Optional marker list used by `run_aegis_full()`.

- ...:

  Additional arguments forwarded to method runners.

## Value

`run_deconvolution()` returns a list with:  
`seu`, `deconv`, `methods_run`, `methods_skipped`, and `messages`.  
`run_aegis_full()` returns an `aegis` object.  
Method-specific runners return spot-by-celltype numeric matrices.
