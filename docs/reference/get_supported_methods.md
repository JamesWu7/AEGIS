# List supported deconvolution methods and execution modes

Returns the AEGIS method registry used by
[`run_deconvolution()`](https://jameswu7.github.io/AEGIS/reference/run_deconvolution.md).
The registry explicitly separates direct R execution, optional Python
execution, and import-only methods.

## Usage

``` r
get_supported_methods()
```

## Value

A `data.frame` with method capabilities, support mode, dependency type,
and adapter/runner mapping.
