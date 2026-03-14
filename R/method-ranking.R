#' Score deconvolution methods from available audit evidence
#'
#' Builds a transparent method-evidence table from existing AEGIS outputs.
#' The final raw score is the mean of normalized evidence dimensions.
#'
#' @param x An `aegis` object.
#' @param use_prior Logical; include prior benchmark scores.
#' @param prior_scores Optional prior scores. Either:
#'   \itemize{
#'     \item named numeric vector (`names` are methods)
#'     \item data.frame with columns `method` and `prior_score`
#'   }
#'
#' @return Modified `aegis` object with `x$consensus$method_evidence`.
#' @export
score_methods <- function(
    x,
    use_prior = FALSE,
    prior_scores = NULL) {
  if (is_multi_sample_context(x)) {
    sample_results <- iterate_aegis_samples(
      x,
      function(obj) score_methods(obj, use_prior = use_prior, prior_scores = prior_scores)
    )
    by_sample <- lapply(sample_results, function(obj) obj$consensus$method_evidence)

    combined <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]
      if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) return(data.frame())
      tbl$sample_id <- sid
      tbl
    }))

    summary_tbl <- if (nrow(combined) > 0L) {
      combined |>
        dplyr::group_by(.data$method) |>
        dplyr::summarise(
          marker_score = mean(.data$marker_score, na.rm = TRUE),
          spatial_score = mean(.data$spatial_score, na.rm = TRUE),
          agreement_score = mean(.data$agreement_score, na.rm = TRUE),
          stability_score = mean(.data$stability_score, na.rm = TRUE),
          prior_score = mean(.data$prior_score, na.rm = TRUE),
          overall_score_raw = mean(.data$overall_score_raw, na.rm = TRUE),
          .groups = "drop"
        )
    } else {
      data.frame()
    }

    x$consensus$method_evidence <- list(
      by_sample = by_sample,
      table = combined,
      summary = summary_tbl,
      use_prior = isTRUE(use_prior)
    )
    return(x)
  }

  assert_is_aegis(x)
  if (!is.logical(use_prior) || length(use_prior) != 1L || is.na(use_prior)) {
    stop("`use_prior` must be TRUE/FALSE.", call. = FALSE)
  }
  if (is.null(x$deconv) || !is.list(x$deconv) || length(x$deconv) == 0L) {
    stop("`x$deconv` must be a non-empty list.", call. = FALSE)
  }

  methods <- names(x$deconv)
  evidence <- data.frame(
    method = methods,
    marker_score = NA_real_,
    spatial_score = NA_real_,
    agreement_score = NA_real_,
    stability_score = NA_real_,
    prior_score = NA_real_,
    stringsAsFactors = FALSE
  )

  has_any_evidence <- FALSE

  # Marker evidence: method-level concordance (prefer Spearman; fallback Pearson).
  marker_summary <- x$audit$marker$summary_table %||% x$audit$marker$summary
  if (!is.null(marker_summary) && is.data.frame(marker_summary) && nrow(marker_summary) > 0L &&
    all(c("method") %in% colnames(marker_summary))) {
    marker_col <- if ("mean_spearman" %in% colnames(marker_summary)) {
      "mean_spearman"
    } else if ("mean_pearson" %in% colnames(marker_summary)) {
      "mean_pearson"
    } else {
      NULL
    }
    if (!is.null(marker_col)) {
      evidence$marker_score <- marker_summary[[marker_col]][match(methods, marker_summary$method)]
      has_any_evidence <- TRUE
    }
  }

  # Spatial evidence: prefer mean_smoothness; fallback inverse inconsistency.
  spatial_summary <- x$audit$spatial$summary
  if (!is.null(spatial_summary) && is.data.frame(spatial_summary) && nrow(spatial_summary) > 0L &&
    all(c("method") %in% colnames(spatial_summary))) {
    if ("mean_smoothness" %in% colnames(spatial_summary)) {
      evidence$spatial_score <- spatial_summary$mean_smoothness[match(methods, spatial_summary$method)]
      has_any_evidence <- TRUE
    } else if ("mean_local_inconsistency" %in% colnames(spatial_summary)) {
      vals <- spatial_summary$mean_local_inconsistency[match(methods, spatial_summary$method)]
      evidence$spatial_score <- 1 / (1 + vals)
      has_any_evidence <- TRUE
    }
  }

  # Agreement evidence: average pairwise correlation participation.
  comp_pairwise <- x$consensus$comparison$pairwise_summary
  if (!is.null(comp_pairwise) && is.data.frame(comp_pairwise) && nrow(comp_pairwise) > 0L &&
    all(c("method_1", "method_2", "mean_correlation") %in% colnames(comp_pairwise))) {
    agreement_vals <- vapply(methods, function(m) {
      idx <- comp_pairwise$method_1 == m | comp_pairwise$method_2 == m
      if (!any(idx)) return(NA_real_)
      mean(comp_pairwise$mean_correlation[idx], na.rm = TRUE)
    }, numeric(1))
    evidence$agreement_score <- agreement_vals
    has_any_evidence <- TRUE
  }

  # Stability evidence from basic audit: penalize row-sum deviation.
  basic_summary <- x$audit$basic$summary_table %||% x$audit$basic$summary
  if (!is.null(basic_summary) && is.data.frame(basic_summary) && nrow(basic_summary) > 0L &&
    all(c("method") %in% colnames(basic_summary))) {
    if ("mean_sum_dev" %in% colnames(basic_summary)) {
      vals <- basic_summary$mean_sum_dev[match(methods, basic_summary$method)]
      evidence$stability_score <- 1 / (1 + vals)
      has_any_evidence <- TRUE
    } else if ("mean_dominance" %in% colnames(basic_summary)) {
      vals <- basic_summary$mean_dominance[match(methods, basic_summary$method)]
      evidence$stability_score <- vals
      has_any_evidence <- TRUE
    }
  }

  if (isTRUE(use_prior)) {
    evidence$prior_score <- parse_prior_scores(prior_scores = prior_scores, methods = methods)
  }

  if (!has_any_evidence && !isTRUE(use_prior)) {
    stop(
      "No method-evidence sources found. Run one or more of: audit_basic(), audit_marker(), audit_spatial(), compare_methods().",
      call. = FALSE
    )
  }

  dims <- c("marker_score", "spatial_score", "agreement_score", "stability_score")
  if (isTRUE(use_prior)) dims <- c(dims, "prior_score")

  norm_cols <- lapply(dims, function(col) normalize_evidence_to_unit(evidence[[col]]))
  norm_mat <- as.data.frame(norm_cols, stringsAsFactors = FALSE)
  names(norm_mat) <- dims

  evidence$overall_score_raw <- rowMeans(norm_mat, na.rm = TRUE)
  evidence$overall_score_raw[rowSums(is.finite(as.matrix(norm_mat))) == 0L] <- NA_real_

  if (all(!is.finite(evidence$overall_score_raw))) {
    stop("Method evidence exists but could not produce finite scores. Check upstream audit outputs.", call. = FALSE)
  }

  x$consensus$method_evidence <- evidence
  x
}

#' Rank methods from evidence dimensions
#'
#' Converts each evidence dimension into ranks and aggregates them with either
#' Robust Rank Aggregation (RRA) or simple mean-rank.
#'
#' @param x An `aegis` object.
#' @param method Aggregation method: `rra` or `mean_rank`.
#' @param use_prior Logical; include prior rank dimension when available.
#'
#' @return Modified `aegis` object with `x$consensus$method_ranking`.
#' @export
rank_methods <- function(
    x,
    method = c("rra", "mean_rank"),
    use_prior = FALSE) {
  method <- match.arg(method)

  if (is_multi_sample_context(x)) {
    sample_results <- iterate_aegis_samples(
      x,
      function(obj) rank_methods(obj, method = method, use_prior = use_prior)
    )
    by_sample <- lapply(sample_results, function(obj) get_method_ranking_table(obj))

    combined <- dplyr::bind_rows(lapply(names(by_sample), function(sid) {
      tbl <- by_sample[[sid]]
      if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0L) return(data.frame())
      tbl$sample_id <- sid
      tbl
    }))

    summary_tbl <- if (nrow(combined) > 0L) {
      combined |>
        dplyr::group_by(.data$method) |>
        dplyr::summarise(
          overall_rank = mean(.data$overall_rank, na.rm = TRUE),
          overall_score = mean(.data$overall_score, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::arrange(.data$overall_rank)
    } else {
      data.frame()
    }
    if (nrow(summary_tbl) > 0L) {
      summary_tbl$recommendation <- recommendation_from_rank(summary_tbl$overall_rank)
    }

    x$consensus$method_ranking <- list(
      by_sample = by_sample,
      table = combined,
      summary = summary_tbl,
      requested_method = method
    )
    return(x)
  }

  assert_is_aegis(x)
  if (!is.logical(use_prior) || length(use_prior) != 1L || is.na(use_prior)) {
    stop("`use_prior` must be TRUE/FALSE.", call. = FALSE)
  }

  evidence <- get_method_evidence_table(x)
  if (is.null(evidence) || !is.data.frame(evidence) || nrow(evidence) == 0L) {
    stop("Method evidence not found. Run score_methods() first.", call. = FALSE)
  }

  dims <- c("marker_score", "spatial_score", "agreement_score", "stability_score")
  if (isTRUE(use_prior)) dims <- c(dims, "prior_score")
  dims <- intersect(dims, colnames(evidence))

  usable_dims <- dims[vapply(dims, function(col) sum(is.finite(evidence[[col]])) >= 2L, logical(1))]
  if (length(usable_dims) == 0L) {
    stop("No usable evidence dimensions to rank methods.", call. = FALSE)
  }

  rank_cols <- lapply(usable_dims, function(col) rank(-evidence[[col]], ties.method = "average", na.last = "keep"))
  rank_tbl <- as.data.frame(rank_cols, stringsAsFactors = FALSE)
  names(rank_tbl) <- paste0("rank_", usable_dims)

  out <- cbind(
    data.frame(method = evidence$method, stringsAsFactors = FALSE),
    rank_tbl
  )

  used_method <- method
  if (identical(method, "rra")) {
    rra_ok <- requireNamespace("RobustRankAggreg", quietly = TRUE)
    rra_tbl <- NULL
    if (isTRUE(rra_ok)) {
      rank_lists <- lapply(names(rank_tbl), function(rc) {
        idx <- order(out[[rc]], na.last = NA)
        out$method[idx]
      })
      rank_lists <- rank_lists[lengths(rank_lists) > 0L]
      if (length(rank_lists) > 0L) {
        rra_tbl <- tryCatch(
          RobustRankAggreg::aggregateRanks(glist = rank_lists),
          error = function(e) NULL
        )
      }
    }

    if (!is.null(rra_tbl) && is.data.frame(rra_tbl) && nrow(rra_tbl) > 0L &&
      all(c("Name", "Score") %in% colnames(rra_tbl))) {
      map_idx <- match(out$method, as.character(rra_tbl$Name))
      out$rra_pvalue <- as.numeric(rra_tbl$Score)[map_idx]
      out$overall_rank <- rank(out$rra_pvalue, ties.method = "average", na.last = "keep")
      out$overall_score <- -log10(pmax(out$rra_pvalue, 1e-16))
    } else {
      warning(
        "RRA aggregation unavailable (missing RobustRankAggreg or failed aggregation); falling back to mean_rank.",
        call. = FALSE
      )
      used_method <- "mean_rank_fallback"
      out$overall_rank <- rowMeans(rank_tbl, na.rm = TRUE)
      out$overall_score <- -out$overall_rank
      out$rra_pvalue <- NA_real_
    }
  } else {
    out$overall_rank <- rowMeans(rank_tbl, na.rm = TRUE)
    out$overall_score <- -out$overall_rank
    out$rra_pvalue <- NA_real_
  }

  out <- out[order(out$overall_rank), , drop = FALSE]
  out$recommendation <- recommendation_from_rank(out$overall_rank)
  out$aggregation_requested <- method
  out$aggregation_used <- used_method

  x$consensus$method_ranking <- out
  x
}

#' @keywords internal
normalize_evidence_to_unit <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  vals <- x[ok]
  rng <- range(vals)
  if (!is.finite(rng[[1]]) || !is.finite(rng[[2]])) {
    out[ok] <- NA_real_
    return(out)
  }
  if (abs(rng[[2]] - rng[[1]]) < .Machine$double.eps) {
    out[ok] <- 0.5
    return(out)
  }
  out[ok] <- (vals - rng[[1]]) / (rng[[2]] - rng[[1]])
  out
}

#' @keywords internal
parse_prior_scores <- function(prior_scores, methods) {
  out <- rep(NA_real_, length(methods))
  if (is.null(prior_scores)) {
    warning("`use_prior = TRUE` but `prior_scores` is NULL; prior dimension will be NA.", call. = FALSE)
    return(out)
  }

  if (is.numeric(prior_scores)) {
    nm <- names(prior_scores)
    if (is.null(nm) || any(trimws(nm) == "")) {
      stop("When numeric, `prior_scores` must be a named vector keyed by method.", call. = FALSE)
    }
    out <- as.numeric(prior_scores[match(methods, nm)])
    return(out)
  }

  if (is.data.frame(prior_scores)) {
    if (!all(c("method", "prior_score") %in% colnames(prior_scores))) {
      stop("`prior_scores` data.frame must include columns `method` and `prior_score`.", call. = FALSE)
    }
    out <- as.numeric(prior_scores$prior_score[match(methods, prior_scores$method)])
    return(out)
  }

  stop("`prior_scores` must be NULL, named numeric vector, or data.frame(method, prior_score).", call. = FALSE)
}

#' @keywords internal
recommendation_from_rank <- function(overall_rank) {
  out <- rep("use_with_caution", length(overall_rank))
  ok <- is.finite(overall_rank)
  n <- sum(ok)
  if (n == 0L) return(out)

  ranks <- rank(overall_rank[ok], ties.method = "average")
  preferred_cut <- max(1L, ceiling(n * 0.33))
  acceptable_cut <- max(preferred_cut, ceiling(n * 0.67))

  out_ok <- rep("use_with_caution", n)
  out_ok[ranks <= preferred_cut] <- "preferred"
  out_ok[ranks > preferred_cut & ranks <= acceptable_cut] <- "acceptable"
  out[ok] <- out_ok
  out
}

#' @keywords internal
get_method_evidence_table <- function(x) {
  me <- x$consensus$method_evidence
  if (is.null(me)) return(NULL)
  if (is.data.frame(me)) return(me)
  if (is.list(me)) {
    if (is.data.frame(me$summary) && nrow(me$summary) > 0L) return(me$summary)
    if (is.data.frame(me$table) && nrow(me$table) > 0L) return(me$table)
  }
  NULL
}

#' @keywords internal
get_method_ranking_table <- function(x) {
  mr <- x$consensus$method_ranking
  if (is.null(mr)) return(NULL)
  if (is.data.frame(mr)) return(mr)
  if (is.list(mr)) {
    if (is.data.frame(mr$summary) && nrow(mr$summary) > 0L) return(mr$summary)
    if (is.data.frame(mr$table) && nrow(mr$table) > 0L) return(mr$table)
  }
  NULL
}
