# AEGIS One-Step Deconvolution Workflow

This tutorial demonstrates one-step method dispatch and downstream AEGIS analysis.

## Step 1. Method registry and dispatch

```r
registry <- get_supported_methods()
res <- run_deconvolution(
  seu = seu,
  reference = reference,
  methods = c("SPOTlight", "CARD", "RCTD"),
  strict = FALSE
)
obj <- run_aegis(res$seu, deconv = res$deconv, markers = markers)
obj <- score_methods(obj)
obj <- rank_methods(obj, method = "mean_rank")
obj <- compute_consensus(obj, strategy = "weighted", top_n = 2)
```

## Step 2. Visualization

```r
plot_compare(obj, type = "heatmap")
```

![heatmap](AEGIS-one-step-deconvolution_files/figure-html/one-step-heatmap.png)

```r
plot_compare(obj, type = "consensus_map")
```

![consensus](AEGIS-one-step-deconvolution_files/figure-html/one-step-consensus.png)

```r
plot_compare(obj, type = "ranking")
```

![ranking](AEGIS-one-step-deconvolution_files/figure-html/one-step-ranking.png)

```r
plot_compare(obj, type = "disagreement_map")
```

![disagreement](AEGIS-one-step-deconvolution_files/figure-html/one-step-disagreement.png)

```r
plot_compare(obj, type = "confidence_map")
```

![confidence](AEGIS-one-step-deconvolution_files/figure-html/one-step-confidence.png)

## Practical notes

- Use R-native runnable methods first.
- If a method is not executable in your environment, use `read_*()` import adapters.
- Keep `plot_compare()` as the primary visualization entry.
