# AEGIS Overview

## Purpose

AEGIS is designed to audit and compare spatial deconvolution outputs on
Seurat spatial objects.

Current MVP workflows:

1.  Simulate deconvolution outputs for demos and method development.
2.  Import real exported results through method adapters
    (RCTD/SPOTlight/cell2location/CARD/DestVI/STdeconvolve and others).
3.  Run a unified audit/comparison/consensus pipeline on the
    standardized matrices.
4.  Run single-sample or multi-sample analysis through one entry point
    ([`run_aegis()`](https://jameswu7.github.io/AEGIS/reference/run_aegis.md)).

## AEGIS Object Structure

[`as_aegis()`](https://jameswu7.github.io/AEGIS/reference/as_aegis.md)
creates a lightweight S3 object:

- `seu`: Seurat object
- `deconv`: named list of spot-by-celltype matrices
- `markers`: optional marker list
- `audit`: audit results (`basic`, `marker`, `spatial`)
- `consensus`: comparison + consensus outputs
- `meta`: user metadata list

## Simulated Workflow

``` r
data("aegis_example", package = "AEGIS")

deconv <- simulate_deconv_results(aegis_example, seed = 123)
#> Loading required namespace: SeuratObject
obj <- run_aegis(
  aegis_example,
  deconv = deconv,
  markers = aegis_default_markers()
)

knitr::kable(obj$audit$basic$summary)
```

| method        | n_spots | n_celltypes | zero_fraction | near_zero_fraction | mean_dominance | mean_entropy | mean_n_detected_types | mean_sum_dev |
|:--------------|--------:|------------:|--------------:|-------------------:|---------------:|-------------:|----------------------:|-------------:|
| RCTD          |    1200 |           7 |     0.1100000 |          0.2673810 |      0.3741816 |     1.547001 |              5.128333 |            0 |
| SPOTlight     |    1200 |           7 |     0.0392857 |          0.1688095 |      0.3092454 |     1.702123 |              5.818333 |            0 |
| cell2location |    1200 |           7 |     0.0913095 |          0.2263095 |      0.3403734 |     1.615380 |              5.415833 |            0 |

## Real Import Workflow

AEGIS imports exported tables and does not run external deconvolution
backends.

``` r
spots <- colnames(aegis_example)[1:6]

tmp_rctd <- tempfile(fileext = ".csv")
utils::write.csv(
  data.frame(
    barcode = spots,
    B_cell = c(0.5, 0.3, 0.2, 0.6, 0.1, 0.4),
    T_cell = c(0.3, 0.4, 0.4, 0.2, 0.7, 0.4),
    Myeloid = c(0.2, 0.3, 0.4, 0.2, 0.2, 0.2),
    check.names = FALSE
  ),
  tmp_rctd,
  row.names = FALSE
)

rctd <- read_rctd(tmp_rctd)
obj_real <- run_aegis(
  aegis_example[, spots],
  deconv = list(RCTD = rctd),
  do_marker = FALSE,
  do_spatial = FALSE,
  do_compare = FALSE,
  do_consensus = FALSE
)
#> Warning: Not validating Centroids objects
#> Not validating Centroids objects
#> Warning: Not validating FOV objects
#> Not validating FOV objects
#> Not validating FOV objects
#> Not validating FOV objects
#> Not validating FOV objects
#> Not validating FOV objects
#> Warning: Not validating Seurat objects
knitr::kable(obj_real$audit$basic$summary)
```

| method | n_spots | n_celltypes | zero_fraction | near_zero_fraction | mean_dominance | mean_entropy | mean_n_detected_types | mean_sum_dev |
|:-------|--------:|------------:|--------------:|-------------------:|---------------:|-------------:|----------------------:|-------------:|
| RCTD   |       6 |           3 |             0 |                  0 |            0.5 |    0.9967471 |                     3 |            0 |

## Multi-sample Workflow

``` r
spots_all <- colnames(aegis_example)
n_half <- floor(length(spots_all) / 2)

seu_list <- list(
  sample1 = suppressWarnings(aegis_example[, spots_all[seq_len(n_half)]]),
  sample2 = suppressWarnings(aegis_example[, spots_all[seq.int(n_half + 1L, length(spots_all))]])
)

deconv_nested <- list(
  sample1 = simulate_deconv_results(seu_list$sample1, methods = c("RCTD", "SPOTlight"), seed = 11),
  sample2 = simulate_deconv_results(seu_list$sample2, methods = c("RCTD", "SPOTlight"), seed = 12)
)

obj_multi <- run_aegis(seu_list, deconv = deconv_nested, markers = aegis_default_markers())
knitr::kable(head(summarize_by_sample(obj_multi)))
```

| sample_id | n_spots | method    | methods_available | mean_dominance | mean_entropy | mean_local_inconsistency | mean_spot_agreement | mean_consensus_confidence |
|:----------|--------:|:----------|:------------------|---------------:|-------------:|-------------------------:|--------------------:|--------------------------:|
| sample1   |     600 | RCTD      | RCTD;SPOTlight    |      0.3771202 |     1.545049 |                0.0940288 |           0.9738783 |                 0.9634940 |
| sample1   |     600 | SPOTlight | RCTD;SPOTlight    |      0.3098760 |     1.703894 |                0.0719510 |           0.9738783 |                 0.9634940 |
| sample2   |     600 | RCTD      | RCTD;SPOTlight    |      0.3774935 |     1.550962 |                0.0958053 |           0.9741903 |                 0.9639251 |
| sample2   |     600 | SPOTlight | RCTD;SPOTlight    |      0.3085310 |     1.700376 |                0.0739535 |           0.9741903 |                 0.9639251 |

## Next Steps

See the demo and complete tutorials for full end-to-end examples,
including multi-sample reporting.
