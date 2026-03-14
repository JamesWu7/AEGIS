#' Compute consensus deconvolution profile
#'
#' Aggregates deconvolution outputs across methods and computes disagreement and
#' confidence summaries.
#'
#' @param x An `aegis` object.
#' @param strategy Aggregation strategy: `mean`, `weighted`, or `trimmed_mean`.
#' @param top_n Optional number of top-ranked methods to include.
#'
#' @return Modified `aegis` object with `x$consensus$result` populated.
#' @export
compute_consensus <- function(
    x,
    strategy = c("mean", "weighted", "trimmed_mean"),
    top_n = NULL) {
  strategy <- match.arg(strategy)

  if (is_multi_sample_context(x)) {
    sample_results <- iterate_aegis_samples(
      x,
      function(obj) compute_consensus(obj, strategy = strategy, top_n = top_n)
    )
    by_sample <- lapply(sample_results, function(obj) obj$consensus$result)

    spot_confidence <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]$spot_confidence %||% data.frame()
      if (nrow(tbl) == 0L) return(tbl)
      tbl$sample_id <- sid
      tbl
    }))
    celltype_stability <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]$celltype_stability %||% data.frame()
      if (nrow(tbl) == 0L) return(tbl)
      tbl$sample_id <- sid
      tbl
    }))
    shared_celltypes <- lapply(by_sample, function(z) z$shared_celltypes %||% character())
    methods_used <- lapply(by_sample, function(z) z$methods_used %||% character())

    x$consensus$result <- list(
      by_sample = by_sample,
      spot_confidence = spot_confidence,
      celltype_stability = celltype_stability,
      shared_celltypes = shared_celltypes,
      methods_used = methods_used,
      strategy_used = strategy,
      top_n = top_n
    )
    return(x)
  }

  assert_is_aegis(x)
  if (is.null(x$deconv) || !is.list(x$deconv) || length(x$deconv) == 0L) {
    stop("`x$deconv` must be a non-empty list.", call. = FALSE)
  }

  methods <- names(x$deconv)
  spots <- colnames(x$seu)
  if (is.null(spots) || length(spots) == 0L) {
    stop("Seurat object has no spot names.", call. = FALSE)
  }
  for (method in methods) {
    mat <- x$deconv[[method]]
    if (!is.matrix(mat) || !is.numeric(mat) || nrow(mat) == 0L || ncol(mat) == 0L) {
      stop(sprintf("Method '%s' must be a non-empty numeric matrix.", method), call. = FALSE)
    }
    if (!identical(rownames(mat), spots)) {
      stop(sprintf("Method '%s' rownames are not aligned to Seurat spots.", method), call. = FALSE)
    }
  }

  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n <= 0 || top_n %% 1 != 0) {
      stop("`top_n` must be NULL or a positive integer.", call. = FALSE)
    }
    top_n <- as.integer(top_n)
    if (top_n > length(methods)) {
      stop(sprintf("`top_n` must be <= number of available methods (%d).", length(methods)), call. = FALSE)
    }

    ranking_tbl <- get_method_ranking_table(x)
    if (is.null(ranking_tbl) || !("method" %in% colnames(ranking_tbl)) || !("overall_rank" %in% colnames(ranking_tbl))) {
      stop("`top_n` requires method ranking. Run rank_methods() first.", call. = FALSE)
    }
    ranking_tbl <- ranking_tbl[ranking_tbl$method %in% methods, , drop = FALSE]
    ranking_tbl <- ranking_tbl[order(ranking_tbl$overall_rank), , drop = FALSE]
    if (nrow(ranking_tbl) < top_n) {
      stop("Not enough ranked methods to satisfy `top_n`.", call. = FALSE)
    }
    methods <- ranking_tbl$method[seq_len(top_n)]
  }

  shared_celltypes <- Reduce(intersect, lapply(x$deconv[methods], colnames))
  if (length(shared_celltypes) == 0L) {
    stop("No shared cell types found across methods for consensus.", call. = FALSE)
  }
  shared_celltypes <- sort(shared_celltypes)

  values_by_method <- lapply(methods, function(method) {
    x$deconv[[method]][spots, shared_celltypes, drop = FALSE]
  })
  names(values_by_method) <- methods
  arr <- simplify2array(values_by_method)

  weights <- rep(1 / length(methods), length(methods))
  names(weights) <- methods
  if (identical(strategy, "weighted")) {
    ranking_tbl <- get_method_ranking_table(x)
    evidence_tbl <- get_method_evidence_table(x)
    if (is.null(ranking_tbl) && is.null(evidence_tbl)) {
      x <- score_methods(x, use_prior = FALSE)
      x <- rank_methods(x, method = "mean_rank", use_prior = FALSE)
      ranking_tbl <- get_method_ranking_table(x)
      evidence_tbl <- get_method_evidence_table(x)
    }

    score_vec <- NULL
    if (!is.null(ranking_tbl) && all(c("method", "overall_score") %in% colnames(ranking_tbl))) {
      score_vec <- ranking_tbl$overall_score[match(methods, ranking_tbl$method)]
    }
    if (is.null(score_vec) || all(!is.finite(score_vec))) {
      if (!is.null(evidence_tbl) && all(c("method", "overall_score_raw") %in% colnames(evidence_tbl))) {
        score_vec <- evidence_tbl$overall_score_raw[match(methods, evidence_tbl$method)]
      }
    }
    if (is.null(score_vec)) {
      score_vec <- rep(1, length(methods))
    }

    score_vec[!is.finite(score_vec)] <- NA_real_
    if (all(is.na(score_vec))) {
      score_vec <- rep(1, length(methods))
    } else {
      min_ok <- min(score_vec, na.rm = TRUE)
      score_vec <- score_vec - min_ok
      score_vec[score_vec <= 0 | !is.finite(score_vec)] <- NA_real_
      if (all(is.na(score_vec))) score_vec <- rep(1, length(methods))
    }
    score_vec[is.na(score_vec)] <- stats::median(score_vec, na.rm = TRUE)
    if (!all(is.finite(score_vec)) || sum(score_vec) <= 0) {
      score_vec <- rep(1, length(methods))
    }
    weights <- score_vec / sum(score_vec)
    names(weights) <- methods
  }

  if (identical(strategy, "trimmed_mean")) {
    if (length(methods) < 3L) {
      warning("trimmed_mean requires at least 3 methods; falling back to mean.", call. = FALSE)
      strategy <- "mean"
    }
  }

  if (identical(strategy, "mean")) {
    consensus <- apply(arr, c(1, 2), mean)
  } else if (identical(strategy, "weighted")) {
    consensus <- apply(arr, c(1, 2), function(v) stats::weighted.mean(v, w = weights, na.rm = TRUE))
  } else {
    trim <- if (length(methods) >= 5L) 0.2 else 0.1
    consensus <- apply(arr, c(1, 2), function(v) mean(v, trim = trim, na.rm = TRUE))
  }

  disagreement <- apply(arr, c(1, 2), function(v) {
    if (!is.finite(sum(v))) return(NA_real_)
    if (length(v) < 2L) return(0)
    if (identical(strategy, "weighted")) {
      weighted_sd(v, weights)
    } else {
      stats::sd(v)
    }
  })

  rownames(consensus) <- spots
  colnames(consensus) <- shared_celltypes
  rownames(disagreement) <- spots
  colnames(disagreement) <- shared_celltypes

  consensus[!is.finite(consensus)] <- NA_real_
  disagreement[!is.finite(disagreement)] <- NA_real_

  rs <- rowSums(consensus)
  rs[!is.finite(rs) | rs <= 0] <- 1
  consensus <- consensus / rs

  mean_sd <- rowMeans(disagreement, na.rm = TRUE)
  mean_sd[!is.finite(mean_sd)] <- NA_real_
  spot_confidence <- data.frame(
    spot = spots,
    mean_disagreement = mean_sd,
    confidence = 1 / (1 + mean_sd),
    stringsAsFactors = FALSE
  )

  ct_sd <- colMeans(disagreement, na.rm = TRUE)
  ct_sd[!is.finite(ct_sd)] <- NA_real_
  celltype_stability <- data.frame(
    celltype = shared_celltypes,
    mean_disagreement = ct_sd,
    stability = 1 / (1 + ct_sd),
    stringsAsFactors = FALSE
  )

  x$consensus$result <- list(
    consensus_matrix = consensus,
    method_disagreement = disagreement,
    spot_confidence = spot_confidence,
    celltype_stability = celltype_stability,
    shared_celltypes = shared_celltypes,
    methods = methods,
    methods_used = methods,
    strategy_used = strategy,
    top_n = top_n,
    method_weights = data.frame(method = methods, weight = as.numeric(weights), stringsAsFactors = FALSE)
  )

  x
}

#' @keywords internal
weighted_sd <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  x <- x[ok]
  w <- w[ok]
  if (length(x) <= 1L) return(0)
  w <- w / sum(w)
  mu <- sum(w * x)
  sqrt(sum(w * (x - mu)^2))
}
