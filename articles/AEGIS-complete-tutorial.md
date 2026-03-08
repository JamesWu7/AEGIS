# AEGIS Complete Tutorial

## 1. Load example Seurat object

``` r
data("aegis_example", package = "AEGIS")
seu <- aegis_example
```

## 2. Simulate method outputs

``` r
deconv <- simulate_deconv_results(
  seu,
  methods = c("RCTD", "SPOTlight", "cell2location"),
  seed = 2026
)
```

## 3. Construct aegis object

``` r
obj <- as_aegis(seu, deconv)
obj
#> <aegis>
#>   Spots: 1200
#>   Genes: 36601
#>   Methods: 3 (RCTD, SPOTlight, cell2location)
#>   Marker sets: none
#>   Audits: none
```

## 4. Run basic audit

``` r
obj <- audit_basic(obj)
knitr::kable(obj$audit$basic$summary)
```

| method        | n_spots | n_celltypes | zero_fraction | near_zero_fraction | mean_dominance | mean_entropy | mean_n_detected_types | mean_sum_dev |
|:--------------|--------:|------------:|--------------:|-------------------:|---------------:|-------------:|----------------------:|-------------:|
| RCTD          |    1200 |           7 |     0.1015476 |          0.2575000 |      0.3695554 |     1.557445 |              5.197500 |            0 |
| SPOTlight     |    1200 |           7 |     0.0410714 |          0.1664286 |      0.3073309 |     1.702981 |              5.835000 |            0 |
| cell2location |    1200 |           7 |     0.0896429 |          0.2244048 |      0.3407826 |     1.617544 |              5.429167 |            0 |

## 5. Inspect spot-level metrics

``` r
head(obj$audit$basic$spot_metrics)
#>                                      spot method      sum_dev dominance
#> AAACAATCTACTAGCA-1...1 AAACAATCTACTAGCA-1   RCTD 1.110223e-16 0.6302157
#> AAACACCAATAACTGC-1...2 AAACACCAATAACTGC-1   RCTD 0.000000e+00 0.3194496
#> AAACATTTCCCGGATT-1...3 AAACATTTCCCGGATT-1   RCTD 0.000000e+00 0.3418494
#> AAACCGGGTAGGTACC-1...4 AAACCGGGTAGGTACC-1   RCTD 0.000000e+00 0.3862453
#> AAACGAAGAACATACC-1...5 AAACGAAGAACATACC-1   RCTD 1.110223e-16 0.4352017
#> AAACGAGACGGTTGAT-1...6 AAACGAGACGGTTGAT-1   RCTD 0.000000e+00 0.2999716
#>                         entropy n_detected_types
#> AAACAATCTACTAGCA-1...1 1.008560                4
#> AAACACCAATAACTGC-1...2 1.556175                5
#> AAACATTTCCCGGATT-1...3 1.659650                5
#> AAACCGGGTAGGTACC-1...4 1.385390                5
#> AAACGAAGAACATACC-1...5 1.585926                6
#> AAACGAGACGGTTGAT-1...6 1.814070                7
```
