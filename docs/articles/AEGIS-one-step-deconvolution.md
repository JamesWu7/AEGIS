# AEGIS One-Step Deconvolution Workflow

This tutorial demonstrates AEGIS one-step deconvolution orchestration and seamless downstream analysis.

## 1) Inspect method support modes

```r
library(AEGIS)
registry <- get_supported_methods()
knitr::kable(registry)
```

AEGIS distinguishes:
- `run_and_import_r`
- `run_and_import_python`
- `import_only`

## 2) Prepare inputs

```r
data("aegis_example", package = "AEGIS")
seu <- aegis_example
markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))
```

## 3) Run one method directly

```r
res_one <- run_deconvolution(
  seu = seu,
  reference = reference,
  methods = "SPOTlight",
  strict = TRUE
)
obj_one <- run_aegis(res_one$seu, deconv = res_one$deconv, markers = markers)
```

## 4) Run multiple methods

```r
res_multi <- run_deconvolution(
  seu = seu,
  reference = reference,
  methods = c("SPOTlight", "CARD", "RCTD"),
  strict = FALSE
)
obj_multi <- run_aegis(res_multi$seu, deconv = res_multi$deconv, markers = markers)
```

## 5) One-step full wrapper

```r
obj_full <- run_aegis_full(
  seu = seu,
  reference = reference,
  methods = c("SPOTlight", "CARD", "RCTD"),
  markers = markers,
  strict = FALSE
)
```

## 6) Optional Python-backed methods

```r
res_py <- run_deconvolution(
  seu = seu,
  reference = reference,
  methods = c("cell2location", "DestVI", "Tangram"),
  use_python = TRUE,
  strict = FALSE
)
```

## 7) Downstream ranking and consensus

```r
deconv_demo <- simulate_deconv_results(seu, methods = c("SPOTlight", "CARD", "RCTD"), seed = 909)
obj <- run_aegis(seu, deconv = deconv_demo, markers = markers)
obj <- score_methods(obj)
obj <- rank_methods(obj, method = "mean_rank")
obj <- compute_consensus(obj, strategy = "weighted", top_n = 2)

plot_compare(obj, type = "heatmap")
plot_compare(obj, type = "consensus_map")
plot_method_ranking(obj)
```

## 8) Practical notes

- Use direct execution only for methods with ready runtime dependencies.
- Use `read_*()` import adapters for import-only methods.
- Keep Python-backed methods optional and environment-checked.
