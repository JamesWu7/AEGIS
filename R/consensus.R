#' Compute consensus deconvolution profile
#'
#' Aggregates deconvolution outputs across methods (mean aggregation for MVP)
#' and computes disagreement/stability summaries.
#'
#' @param x An `aegis` object.
#'
#' @return Modified `aegis` object with `x$consensus$result` populated.
#' @export
compute_consensus <- function(x) {
  assert_is_aegis(x)

  methods <- names(x$deconv)
  if (length(methods) == 0L) {
    stop("No deconvolution methods available in `x$deconv`.", call. = FALSE)
  }

  spots <- colnames(x$seu)
  all_celltypes <- sort(unique(unlist(lapply(x$deconv, colnames))))

  consensus <- matrix(NA_real_, nrow = length(spots), ncol = length(all_celltypes),
    dimnames = list(spots, all_celltypes)
  )
  disagreement <- consensus

  for (ct in all_celltypes) {
    vals_by_method <- lapply(x$deconv, function(mat) {
      if (ct %in% colnames(mat)) {
        mat[spots, ct]
      } else {
        rep(NA_real_, length(spots))
      }
    })
    vals_mat <- do.call(cbind, vals_by_method)

    consensus[, ct] <- rowMeans(vals_mat, na.rm = TRUE)
    disagreement[, ct] <- apply(vals_mat, 1, stats::sd, na.rm = TRUE)
  }

  consensus[!is.finite(consensus)] <- 0
  disagreement[!is.finite(disagreement)] <- NA_real_

  rs <- rowSums(consensus)
  rs[rs <= 0] <- 1
  consensus <- consensus / rs

  mean_sd <- rowMeans(disagreement, na.rm = TRUE)
  spot_confidence <- data.frame(
    spot = spots,
    mean_sd = mean_sd,
    confidence = 1 / (1 + mean_sd),
    stringsAsFactors = FALSE
  )

  celltype_stability <- data.frame(
    celltype = all_celltypes,
    mean_sd = colMeans(disagreement, na.rm = TRUE),
    stability = 1 / (1 + colMeans(disagreement, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )

  x$consensus$result <- list(
    consensus_matrix = consensus,
    method_disagreement = disagreement,
    spot_confidence = spot_confidence,
    celltype_stability = celltype_stability,
    methods = methods
  )

  x
}
