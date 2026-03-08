# AEGIS Overview

## Purpose

AEGIS provides a minimal, reproducible workflow for spatial
deconvolution audit on Seurat objects.

## P2 workflow

``` r
data("aegis_example", package = "AEGIS")

deconv <- simulate_deconv_results(aegis_example, seed = 123)
#> Loading required namespace: SeuratObject
obj <- as_aegis(aegis_example, deconv)
obj <- audit_basic(obj)

knitr::kable(obj$audit$basic$summary)
```

| method        | n_spots | n_celltypes | zero_fraction | near_zero_fraction | mean_dominance | mean_entropy | mean_n_detected_types | mean_sum_dev |
|:--------------|--------:|------------:|--------------:|-------------------:|---------------:|-------------:|----------------------:|-------------:|
| RCTD          |    1200 |           7 |     0.1100000 |          0.2673810 |      0.3741816 |     1.547001 |              5.128333 |            0 |
| SPOTlight     |    1200 |           7 |     0.0392857 |          0.1688095 |      0.3092454 |     1.702123 |              5.818333 |            0 |
| cell2location |    1200 |           7 |     0.0913095 |          0.2263095 |      0.3403734 |     1.615380 |              5.415833 |            0 |
