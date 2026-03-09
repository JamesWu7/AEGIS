#' Load a set of 10x spatial datasets
#'
#' Loads multiple spatial datasets into a named list of Seurat objects and
#' attaches sample-aware metadata.
#'
#' @param paths Character vector of sample directories.
#' @param sample_ids Optional sample identifiers.
#' @param section_ids Optional section identifiers.
#' @param strict If `TRUE`, ambiguous ID inference errors clearly.
#'
#' @return Named list of Seurat objects.
#' @export
load_10x_spatial_set <- function(paths, sample_ids = NULL, section_ids = NULL, strict = TRUE) {
  if (!is.character(paths) || length(paths) == 0L || anyNA(paths) || any(trimws(paths) == "")) {
    stop("`paths` must be a non-empty character vector of directories.", call. = FALSE)
  }

  paths <- normalizePath(paths, winslash = "/", mustWork = FALSE)
  missing_paths <- paths[!dir.exists(paths)]
  if (length(missing_paths) > 0L) {
    stop(
      sprintf("Some sample paths do not exist: %s", paste(utils::head(missing_paths, 5L), collapse = ", ")),
      call. = FALSE
    )
  }

  inferred_ids <- names(paths)
  if (is.null(sample_ids)) {
    if (!is.null(inferred_ids) && all(nzchar(inferred_ids))) {
      sample_ids <- inferred_ids
    } else {
      sample_ids <- basename(paths)
    }
  }
  if (!is.character(sample_ids) || length(sample_ids) != length(paths)) {
    stop("`sample_ids` must be NULL or a character vector with the same length as `paths`.", call. = FALSE)
  }
  if (anyNA(sample_ids) || any(trimws(sample_ids) == "")) {
    stop("`sample_ids` contains missing or empty names.", call. = FALSE)
  }
  if (anyDuplicated(sample_ids)) {
    stop("`sample_ids` must be unique.", call. = FALSE)
  }

  if (!is.null(section_ids)) {
    if (!is.character(section_ids) || length(section_ids) != length(paths)) {
      stop("`section_ids` must be NULL or a character vector with the same length as `paths`.", call. = FALSE)
    }
  } else {
    section_ids <- sample_ids
  }

  seu_list <- lapply(seq_along(paths), function(i) {
    p <- paths[[i]]
    sid <- sample_ids[[i]]
    sec <- section_ids[[i]]

    seu <- tryCatch({
      # Backward-compatible path for MVP authoritative file naming.
      if (file.exists(file.path(p, "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5"))) {
        load_10x_lymphnode(data_dir = p)
      } else {
        Seurat::Load10X_Spatial(data.dir = p)
      }
    }, error = function(e) {
      stop(
        sprintf("Failed to load sample '%s' from '%s': %s", sid, p, conditionMessage(e)),
        call. = FALSE
      )
    })

    attach_sample_metadata(seu, sample_id = sid, section_id = sec)
  })
  names(seu_list) <- sample_ids

  validate_seurat_list(seu_list, strict = strict)
}

#' Merge a list of spatial Seurat objects
#'
#' @param seu_list Named list of Seurat objects.
#' @param sample_ids Optional sample identifiers.
#' @param add_cell_ids If `TRUE`, prefixes spots with sample IDs before merge.
#' @param preserve_images If `TRUE`, checks image preservation and warns/errors.
#' @param strict If `TRUE`, enforce strict validation.
#'
#' @return A merged Seurat object with `sample_id` metadata.
#' @export
merge_spatial_seurat_list <- function(
    seu_list,
    sample_ids = NULL,
    add_cell_ids = TRUE,
    preserve_images = TRUE,
    strict = TRUE) {
  seu_list <- validate_seurat_list(seu_list, sample_ids = sample_ids, strict = strict)
  sample_ids <- names(seu_list)

  seu_list <- lapply(seq_along(seu_list), function(i) {
    sid <- sample_ids[[i]]
    attach_sample_metadata(seu_list[[i]], sample_id = sid, section_id = sid)
  })
  names(seu_list) <- sample_ids

  if (!isTRUE(add_cell_ids)) {
    all_cells <- unlist(lapply(seu_list, colnames), use.names = FALSE)
    if (anyDuplicated(all_cells)) {
      stop("Spot names are not unique across samples. Set `add_cell_ids = TRUE`.", call. = FALSE)
    }
  }

  merged <- if (length(seu_list) == 1L) {
    seu_list[[1L]]
  } else {
    merge(
      x = seu_list[[1L]],
      y = seu_list[-1L],
      add.cell.ids = if (isTRUE(add_cell_ids)) sample_ids else NULL,
      merge.data = TRUE
    )
  }

  if (!("sample_id" %in% colnames(merged@meta.data))) {
    if (isTRUE(add_cell_ids)) {
      prefix <- vapply(strsplit(colnames(merged), "_", fixed = TRUE), `[[`, character(1), 1L)
      if (!all(prefix %in% sample_ids) && isTRUE(strict)) {
        stop("Could not recover `sample_id` from merged spot names.", call. = FALSE)
      }
      merged$sample_id <- prefix
    } else {
      merged$sample_id <- merged$orig.ident %||% "sample1"
    }
  }

  if (!("section_id" %in% colnames(merged@meta.data))) {
    merged$section_id <- merged$sample_id
  }
  if (!("original_spot" %in% colnames(merged@meta.data))) {
    merged$original_spot <- colnames(merged)
  }

  if (isTRUE(preserve_images)) {
    n_images <- length(merged@images)
    if (n_images == 0L) {
      msg <- "Merged object has no spatial images. Spatial plotting may be unavailable."
      if (isTRUE(strict)) {
        stop(msg, call. = FALSE)
      } else {
        warning(msg, call. = FALSE)
      }
    }
  }

  merged
}

#' Construct a multi-sample AEGIS object
#'
#' @param seu_list Named list of Seurat spatial objects.
#' @param deconv Nested named list of deconvolution matrices by sample.
#' @param markers Optional marker list.
#' @param meta Optional metadata list.
#' @param strict If `TRUE`, enforce strict validation.
#'
#' @return An object of class `aegis_multi`.
#' @export
as_aegis_multi <- function(seu_list, deconv, markers = NULL, meta = NULL, strict = TRUE) {
  seu_list <- validate_seurat_list(seu_list, strict = strict)
  deconv <- validate_nested_deconv_list(deconv, seu_list = seu_list, strict = strict)

  if (!is.null(markers)) {
    if (!is.list(markers)) {
      stop("`markers` must be NULL or a named list.", call. = FALSE)
    }
    marker_names <- names(markers)
    if (length(markers) > 0L && (is.null(marker_names) || any(trimws(marker_names) == ""))) {
      stop("`markers`, when supplied, must be a named list.", call. = FALSE)
    }
  }

  if (is.null(meta)) {
    meta <- list()
  } else if (!is.list(meta)) {
    meta <- tryCatch(as.list(meta), error = function(e) {
      stop("`meta` must be NULL or coercible to list.", call. = FALSE)
    })
  }

  out <- list(
    seu = NULL,
    seu_list = seu_list,
    deconv = deconv,
    markers = markers,
    audit = list(),
    consensus = list(),
    meta = utils::modifyList(
      list(
        mode = "multi",
        sample_ids = names(seu_list),
        n_samples = length(seu_list)
      ),
      meta
    )
  )
  class(out) <- "aegis_multi"
  out
}

#' Split an AEGIS object by sample
#'
#' @param x An `aegis` or `aegis_multi` object.
#'
#' @return Named list of per-sample `aegis` objects.
#' @export
split_aegis_by_sample <- function(x) {
  if (inherits(x, "aegis_multi")) {
    sample_ids <- names(x$seu_list)
    out <- lapply(sample_ids, function(sid) {
      audit_list <- list()
      if (is.list(x$audit) && length(x$audit) > 0L) {
        for (nm in names(x$audit)) {
          if (is.list(x$audit[[nm]]) && "by_sample" %in% names(x$audit[[nm]]) && sid %in% names(x$audit[[nm]]$by_sample)) {
            audit_list[[nm]] <- x$audit[[nm]]$by_sample[[sid]]
          }
        }
      }
      consensus_list <- list()
      if (is.list(x$consensus) && length(x$consensus) > 0L) {
        for (nm in names(x$consensus)) {
          if (is.list(x$consensus[[nm]]) && "by_sample" %in% names(x$consensus[[nm]]) && sid %in% names(x$consensus[[nm]]$by_sample)) {
            consensus_list[[nm]] <- x$consensus[[nm]]$by_sample[[sid]]
          }
        }
      }

      obj <- list(
        seu = x$seu_list[[sid]],
        deconv = x$deconv[[sid]],
        markers = x$markers,
        audit = audit_list,
        consensus = consensus_list,
        meta = utils::modifyList(x$meta %||% list(), list(sample_id = sid, parent_mode = "multi"))
      )
      class(obj) <- "aegis"
      obj
    })
    names(out) <- sample_ids
    return(out)
  }

  assert_is_aegis(x)
  sample_map <- map_spots_to_sample(x$seu)
  sample_ids <- unique(sample_map$sample_id)

  out <- lapply(sample_ids, function(sid) {
    cells <- sample_map$spot[sample_map$sample_id == sid]
    seu_sub <- suppressWarnings(x$seu[, cells])
    deconv_sub <- lapply(names(x$deconv), function(method) {
      mat <- x$deconv[[method]]
      missing <- setdiff(cells, rownames(mat))
      if (length(missing) > 0L) {
        stop(
          sprintf("Method '%s' is missing %d spots for sample '%s'.", method, length(missing), sid),
          call. = FALSE
        )
      }
      mat[cells, , drop = FALSE]
    })
    names(deconv_sub) <- names(x$deconv)

    audit_list <- list()
    if (is.list(x$audit) && length(x$audit) > 0L) {
      for (nm in names(x$audit)) {
        if (is.list(x$audit[[nm]]) && "by_sample" %in% names(x$audit[[nm]]) && sid %in% names(x$audit[[nm]]$by_sample)) {
          audit_list[[nm]] <- x$audit[[nm]]$by_sample[[sid]]
        }
      }
    }
    consensus_list <- list()
    if (is.list(x$consensus) && length(x$consensus) > 0L) {
      for (nm in names(x$consensus)) {
        if (is.list(x$consensus[[nm]]) && "by_sample" %in% names(x$consensus[[nm]]) && sid %in% names(x$consensus[[nm]]$by_sample)) {
          consensus_list[[nm]] <- x$consensus[[nm]]$by_sample[[sid]]
        }
      }
    }

    obj <- list(
      seu = seu_sub,
      deconv = deconv_sub,
      markers = x$markers,
      audit = audit_list,
      consensus = consensus_list,
      meta = utils::modifyList(x$meta %||% list(), list(sample_id = sid, parent_mode = "single_or_merged"))
    )
    class(obj) <- "aegis"
    obj
  })
  names(out) <- sample_ids
  out
}

#' Summarize AEGIS results by sample
#'
#' @param x An `aegis` or `aegis_multi` object.
#'
#' @return A sample-aware summary data.frame.
#' @export
summarize_by_sample <- function(x) {
  sample_objs <- split_aegis_by_sample(x)

  rows <- lapply(names(sample_objs), function(sid) {
    obj <- sample_objs[[sid]]

    if (is.null(obj$audit$basic)) {
      obj <- audit_basic(obj)
    }

    basic <- obj$audit$basic$summary %||% data.frame()
    methods <- names(obj$deconv)
    methods_label <- paste(methods, collapse = ";")

    spatial_tbl <- obj$audit$spatial$summary %||% data.frame()
    comp_tbl <- obj$consensus$comparison$spot_agreement %||% data.frame()
    conf_tbl <- obj$consensus$result$spot_confidence %||% data.frame()

    do.call(rbind, lapply(methods, function(m) {
      basic_m <- basic[basic$method == m, , drop = FALSE]
      spatial_m <- spatial_tbl[spatial_tbl$method == m, , drop = FALSE]
      data.frame(
        sample_id = sid,
        n_spots = ncol(obj$seu),
        method = m,
        methods_available = methods_label,
        mean_dominance = if (nrow(basic_m) == 0L) NA_real_ else basic_m$mean_dominance[[1L]],
        mean_entropy = if (nrow(basic_m) == 0L) NA_real_ else basic_m$mean_entropy[[1L]],
        mean_local_inconsistency = if (nrow(spatial_m) == 0L) NA_real_ else spatial_m$mean_local_inconsistency[[1L]],
        mean_spot_agreement = if (nrow(comp_tbl) == 0L) NA_real_ else mean(comp_tbl$agreement, na.rm = TRUE),
        mean_consensus_confidence = if (nrow(conf_tbl) == 0L) NA_real_ else mean(conf_tbl$confidence, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  })

  dplyr::bind_rows(rows)
}

#' Render one report per sample
#'
#' @param x An `aegis` or `aegis_multi` object.
#' @param output_dir Output directory.
#' @param overwrite Whether to overwrite existing files.
#'
#' @return List with report paths and batch summary path.
#' @export
render_report_batch <- function(x, output_dir = "reports", overwrite = FALSE) {
  if (!is.character(output_dir) || length(output_dir) != 1L || trimws(output_dir) == "") {
    stop("`output_dir` must be a single non-empty path.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop(sprintf("Could not create output directory: %s", output_dir), call. = FALSE)
  }

  sample_objs <- split_aegis_by_sample(x)
  report_paths <- vapply(names(sample_objs), function(sid) {
    out <- file.path(output_dir, sprintf("%s_aegis_report.html", sid))
    if (file.exists(out) && !isTRUE(overwrite)) {
      stop(sprintf("Output already exists (set overwrite = TRUE): %s", out), call. = FALSE)
    }
    render_report(sample_objs[[sid]], output_file = out)
  }, character(1))

  summary_tbl <- summarize_by_sample(x)
  summary_file <- file.path(output_dir, "aegis_batch_summary.csv")
  if (file.exists(summary_file) && !isTRUE(overwrite)) {
    stop(sprintf("Batch summary already exists (set overwrite = TRUE): %s", summary_file), call. = FALSE)
  }
  utils::write.csv(summary_tbl, summary_file, row.names = FALSE)

  list(reports = report_paths, summary_file = summary_file)
}

#' @keywords internal
validate_seurat_list <- function(seu_list, sample_ids = NULL, strict = TRUE) {
  if (!is.list(seu_list) || length(seu_list) == 0L) {
    stop("`seu_list` must be a non-empty list of Seurat objects.", call. = FALSE)
  }
  if (!is.null(sample_ids)) {
    if (!is.character(sample_ids) || length(sample_ids) != length(seu_list)) {
      stop("`sample_ids` must be NULL or a character vector matching `seu_list` length.", call. = FALSE)
    }
    names(seu_list) <- sample_ids
  }

  ids <- names(seu_list)
  if (is.null(ids) || anyNA(ids) || any(trimws(ids) == "")) {
    if (isTRUE(strict)) {
      stop("`seu_list` must be named with non-empty sample IDs.", call. = FALSE)
    }
    ids <- paste0("sample", seq_along(seu_list))
    names(seu_list) <- ids
  }
  if (anyDuplicated(ids)) {
    stop("`seu_list` sample names must be unique.", call. = FALSE)
  }

  bad <- names(seu_list)[!vapply(seu_list, inherits, logical(1), what = "Seurat")]
  if (length(bad) > 0L) {
    stop(sprintf("All elements of `seu_list` must be Seurat objects (bad: %s).", paste(bad, collapse = ", ")), call. = FALSE)
  }

  seu_list
}

#' @keywords internal
validate_nested_deconv_list <- function(deconv, seu_list, strict = TRUE) {
  seu_list <- validate_seurat_list(seu_list, strict = strict)
  sample_ids <- names(seu_list)

  if (!is.list(deconv) || length(deconv) == 0L) {
    stop("`deconv` must be a non-empty nested list by sample.", call. = FALSE)
  }
  if (is.null(names(deconv)) || any(trimws(names(deconv)) == "")) {
    stop("Top-level `deconv` list must be named by sample IDs.", call. = FALSE)
  }

  missing_samples <- setdiff(sample_ids, names(deconv))
  extra_samples <- setdiff(names(deconv), sample_ids)
  if (length(missing_samples) > 0L) {
    stop(sprintf("`deconv` is missing samples: %s", paste(missing_samples, collapse = ", ")), call. = FALSE)
  }
  if (length(extra_samples) > 0L && isTRUE(strict)) {
    stop(sprintf("`deconv` has unknown samples: %s", paste(extra_samples, collapse = ", ")), call. = FALSE)
  }

  out <- lapply(sample_ids, function(sid) {
    d <- deconv[[sid]]
    if (!is.list(d) || length(d) == 0L) {
      stop(sprintf("Sample '%s': deconv entry must be a non-empty named list of methods.", sid), call. = FALSE)
    }
    validate_deconv_list(d, seurat_spots = colnames(seu_list[[sid]]))
  })
  names(out) <- sample_ids
  out
}

#' @keywords internal
attach_sample_metadata <- function(seu, sample_id, section_id = NULL) {
  assert_is_seurat(seu)
  if (!is.character(sample_id) || length(sample_id) != 1L || is.na(sample_id) || trimws(sample_id) == "") {
    stop("`sample_id` must be a single non-empty string.", call. = FALSE)
  }
  if (is.null(section_id)) {
    section_id <- sample_id
  }

  if (!("original_spot" %in% colnames(seu@meta.data))) {
    seu$original_spot <- colnames(seu)
  }
  seu$sample_id <- sample_id
  seu$section_id <- section_id
  seu
}

#' @keywords internal
get_sample_ids <- function(x) {
  if (inherits(x, "aegis_multi")) {
    ids <- x$meta$sample_ids %||% names(x$seu_list)
    return(as.character(ids))
  }
  if (inherits(x, "aegis")) {
    if ("sample_id" %in% colnames(x$seu@meta.data)) {
      ids <- unique(as.character(x$seu$sample_id))
      ids <- ids[!is.na(ids) & trimws(ids) != ""]
      if (length(ids) == 0L) {
        return("sample1")
      }
      return(ids)
    }
    return("sample1")
  }
  stop("`x` must be an aegis or aegis_multi object.", call. = FALSE)
}

#' @keywords internal
align_nested_deconv_to_samples <- function(deconv, seu_list, strict = TRUE) {
  validate_nested_deconv_list(deconv = deconv, seu_list = seu_list, strict = strict)
}

#' @keywords internal
flatten_nested_deconv_if_needed <- function(deconv, sample_id = NULL) {
  if (is.list(deconv) && length(deconv) > 0L && all(vapply(deconv, is.list, logical(1)))) {
    if (is.null(sample_id)) {
      stop("`sample_id` is required when flattening nested deconvolution lists.", call. = FALSE)
    }
    if (!(sample_id %in% names(deconv))) {
      stop(sprintf("Sample '%s' not found in nested deconv list.", sample_id), call. = FALSE)
    }
    return(deconv[[sample_id]])
  }
  deconv
}

#' @keywords internal
map_spots_to_sample <- function(seu) {
  assert_is_seurat(seu)
  spots <- colnames(seu)
  if (length(spots) == 0L) {
    stop("Seurat object has no spots.", call. = FALSE)
  }
  sample_id <- if ("sample_id" %in% colnames(seu@meta.data)) as.character(seu$sample_id) else rep("sample1", length(spots))
  sample_id[is.na(sample_id) | trimws(sample_id) == ""] <- "sample1"
  data.frame(spot = spots, sample_id = sample_id, stringsAsFactors = FALSE)
}

#' @keywords internal
extract_spatial_coords_multi <- function(x) {
  if (inherits(x, "aegis_multi")) {
    out <- lapply(names(x$seu_list), function(sid) {
      coords <- extract_spatial_coords(x$seu_list[[sid]])
      coords$sample_id <- sid
      coords
    })
    return(dplyr::bind_rows(out))
  }

  assert_is_aegis(x)
  coords <- extract_spatial_coords(x$seu)
  spot_map <- map_spots_to_sample(x$seu)
  dplyr::left_join(coords, spot_map, by = c("spot" = "spot"))
}

#' @keywords internal
iterate_aegis_samples <- function(x, fn, ...) {
  objs <- split_aegis_by_sample(x)
  out <- lapply(objs, function(obj) fn(obj, ...))
  names(out) <- names(objs)
  out
}

#' @keywords internal
coerce_single_to_multi_aegis <- function(x) {
  assert_is_aegis(x)
  sample_objs <- split_aegis_by_sample(x)
  seu_list <- lapply(sample_objs, `[[`, "seu")
  deconv_nested <- lapply(sample_objs, `[[`, "deconv")

  out <- list(
    seu = NULL,
    seu_list = seu_list,
    deconv = deconv_nested,
    markers = x$markers,
    audit = x$audit %||% list(),
    consensus = x$consensus %||% list(),
    meta = utils::modifyList(
      x$meta %||% list(),
      list(mode = "multi_from_single", sample_ids = names(sample_objs), n_samples = length(sample_objs))
    )
  )
  class(out) <- "aegis_multi"
  out
}

#' @keywords internal
is_multi_sample_context <- function(x) {
  if (inherits(x, "aegis_multi")) {
    return(TRUE)
  }
  if (inherits(x, "aegis") && "sample_id" %in% colnames(x$seu@meta.data)) {
    return(length(unique(as.character(x$seu$sample_id))) > 1L)
  }
  FALSE
}
