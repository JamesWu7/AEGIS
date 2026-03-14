# Score deconvolution methods from available audit evidence

Builds a transparent method-evidence table from existing AEGIS outputs.
The final raw score is the mean of normalized evidence dimensions.

## Usage

``` r
score_methods(x, use_prior = FALSE, prior_scores = NULL)
```

## Arguments

- x:

  An `aegis` object.

- use_prior:

  Logical; include prior benchmark scores.

- prior_scores:

  Optional prior scores. Either:

  - named numeric vector (`names` are methods)

  - data.frame with columns `method` and `prior_score`

## Value

Modified `aegis` object with `x$consensus$method_evidence`.
