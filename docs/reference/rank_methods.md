# Rank methods from evidence dimensions

Converts each evidence dimension into ranks and aggregates them with
either Robust Rank Aggregation (RRA) or simple mean-rank.

## Usage

``` r
rank_methods(x, method = c("rra", "mean_rank"), use_prior = FALSE)
```

## Arguments

- x:

  An `aegis` object.

- method:

  Aggregation method: `rra` or `mean_rank`.

- use_prior:

  Logical; include prior rank dimension when available.

## Value

Modified `aegis` object with `x$consensus$method_ranking`.
