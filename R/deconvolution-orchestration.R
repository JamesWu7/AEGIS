#' List supported deconvolution methods and execution modes
#'
#' Returns the AEGIS method registry used by `run_deconvolution()`.
#' The registry explicitly distinguishes:
#' \itemize{
#'   \item `run_and_import_r`: runnable from R when method packages are available
#'   \item `run_and_import_python`: runnable via `reticulate` when Python modules are available
#'   \item `import_only`: supported through `read_*()` adapters only
#'   \item `experimental`: adapter exists but direct execution path is not stable yet
#' }
#'
#' @return A `data.frame` with method capabilities and adapter/runner mapping.
#' @export
get_supported_methods <- function() {
  data.frame(
    method_name = c(
      "RCTD", "SPOTlight", "cell2location", "CARD", "SpatialDWLS",
      "stereoscope", "DestVI", "Tangram", "STdeconvolve", "DSTG", "STRIDE"
    ),
    support_mode = c(
      "run_and_import_r", "run_and_import_r", "run_and_import_python", "run_and_import_r", "import_only",
      "run_and_import_python", "run_and_import_python", "run_and_import_python", "import_only", "import_only", "import_only"
    ),
    can_run_in_r = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    can_run_in_python = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
    requires_reference = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
    requires_spatial_coords = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
    expected_output_type = c(
      "proportion", "proportion", "abundance", "proportion", "proportion",
      "proportion", "abundance", "mapping", "latent", "proportion", "proportion"
    ),
    adapter_reader = c(
      "read_rctd", "read_spotlight", "read_cell2location", "read_card", "read_spatialdwls",
      "read_stereoscope", "read_destvi", "read_tangram", "read_stdeconvolve", "read_dstg", "read_stride"
    ),
    reader_function = c(
      "read_rctd", "read_spotlight", "read_cell2location", "read_card", "read_spatialdwls",
      "read_stereoscope", "read_destvi", "read_tangram", "read_stdeconvolve", "read_dstg", "read_stride"
    ),
    runner_function = c(
      "run_rctd", "run_spotlight", "run_cell2location", "run_card", NA,
      "run_stereoscope", "run_destvi", "run_tangram", NA, NA, NA
    ),
    dependency_type = c(
      "r_package", "r_package", "python_module", "r_package", "import_only",
      "python_module", "python_module", "python_module", "import_only", "import_only", "import_only"
    ),
    backend_dependency = c(
      "spacexr", "SPOTlight", "cell2location,scvi", "CARD", NA,
      "scvi", "scvi", "tangram", NA, NA, NA
    ),
    notes = c(
      "R-native runner; prefers spacexr pipeline when available.",
      "R-native runner; requires SPOTlight package and suitable reference inputs.",
      "Python via reticulate; if unavailable, import exported tables with read_cell2location().",
      "R-native runner; requires CARD package and suitable reference inputs.",
      "Import-only in P9. Use read_spatialdwls().",
      "Python via reticulate (experimental wrapper).",
      "Python via reticulate (experimental wrapper).",
      "Python via reticulate (experimental wrapper).",
      "Import-only in P9. Use read_stdeconvolve().",
      "Import-only in P9. Use read_dstg().",
      "Import-only in P9. Use read_stride()."
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

#' Run selected deconvolution methods in one unified interface
#'
#' Dispatches methods by registry capabilities and returns standardized
#' spot-by-celltype matrices that can directly enter `as_aegis()` / `run_aegis()`.
#'
#' @param seu A Seurat spatial object.
#' @param reference Reference input used by runnable methods.
#' @param methods Character vector of method names.
#' @param sample_id Optional sample identifier recorded in output metadata.
#' @param normalize Logical; row-normalize outputs to proportions when possible.
#' @param strict Logical; if `TRUE`, fail on any requested method that cannot run.
#' @param use_python Logical; allow Python-backed methods (`reticulate`).
#' @param runner_overrides Optional named list of runner functions keyed by method.
#'   Intended for advanced usage/testing.
#' @param runner_args Optional named list of per-method argument lists passed to
#'   the selected runner.
#' @param ... Additional arguments forwarded to method runners.
#'
#' @return A list with `seu`, `deconv`, `methods_run`, `methods_skipped`, and `messages`.
#' @export
run_deconvolution <- function(
    seu,
    reference,
    methods = c("SPOTlight", "RCTD", "CARD"),
    sample_id = NULL,
    normalize = TRUE,
    strict = TRUE,
    use_python = TRUE,
    runner_overrides = NULL,
    runner_args = NULL,
    ...) {
  assert_is_seurat(seu, "seu")

  if (!is.character(methods) || length(methods) == 0L || anyNA(methods) || any(trimws(methods) == "")) {
    stop("`methods` must be a non-empty character vector.", call. = FALSE)
  }
  methods <- unique(as.character(methods))

  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    stop("`normalize` must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE/FALSE.", call. = FALSE)
  }
  if (!is.logical(use_python) || length(use_python) != 1L || is.na(use_python)) {
    stop("`use_python` must be TRUE/FALSE.", call. = FALSE)
  }

  registry <- get_supported_methods()
  resolved <- resolve_method_names(methods, registry, unique_out = TRUE)

  if (!is.null(runner_overrides)) {
    if (!is.list(runner_overrides) || is.null(names(runner_overrides)) || any(trimws(names(runner_overrides)) == "")) {
      stop("`runner_overrides` must be NULL or a named list of functions.", call. = FALSE)
    }
    names(runner_overrides) <- resolve_method_names(names(runner_overrides), registry = registry, unique_out = FALSE)
  }
  if (!is.null(runner_args) && !is.list(runner_args)) {
    stop("`runner_args` must be NULL or a named list.", call. = FALSE)
  }
  if (!is.null(runner_args) && !is.null(names(runner_args))) {
    if (any(trimws(names(runner_args)) == "")) {
      stop("`runner_args` names must be non-empty method names.", call. = FALSE)
    }
    names(runner_args) <- resolve_method_names(names(runner_args), registry = registry, unique_out = FALSE)
  }

  dots <- list(...)
  spots <- colnames(seu)
  if (is.null(spots) || length(spots) == 0L) {
    stop("`seu` must contain spot names in `colnames(seu)`.", call. = FALSE)
  }

  deconv_out <- list()
  skipped <- data.frame(method = character(0), reason = character(0), stringsAsFactors = FALSE)
  messages <- character(0)

  for (method in resolved) {
    reg_row <- registry[registry$method_name == method, , drop = FALSE]
    mode <- reg_row$support_mode[[1L]]

    if (identical(mode, "import_only")) {
      msg <- sprintf(
        "Method '%s' is import_only in AEGIS. Use %s() and pass results via `deconv`.",
        method,
        reg_row$adapter_reader[[1L]]
      )
      if (isTRUE(strict)) stop(msg, call. = FALSE)
      skipped <- rbind(skipped, data.frame(method = method, reason = "import_only", stringsAsFactors = FALSE))
      messages <- c(messages, msg)
      next
    }

    if (isTRUE(reg_row$requires_reference[[1L]])) {
      validate_reference_input(reference = reference, method = method)
    }

    if (identical(mode, "run_and_import_python") && !isTRUE(use_python)) {
      msg <- sprintf("Method '%s' requires Python execution but `use_python = FALSE`.", method)
      if (isTRUE(strict)) stop(msg, call. = FALSE)
      skipped <- rbind(skipped, data.frame(method = method, reason = "python_disabled", stringsAsFactors = FALSE))
      messages <- c(messages, msg)
      next
    }

    runner <- pick_method_runner(
      method = method,
      registry = registry,
      runner_overrides = runner_overrides
    )

    method_args <- if (!is.null(runner_args) && method %in% names(runner_args) && is.list(runner_args[[method]])) {
      runner_args[[method]]
    } else {
      list()
    }

    runner_call <- c(
      list(
        seu = seu,
        reference = reference,
        normalize = normalize,
        strict = strict
      ),
      method_args,
      dots
    )

    mat <- tryCatch(
      do.call(runner, runner_call),
      error = function(e) {
        msg <- sprintf("Method '%s' execution failed: %s", method, conditionMessage(e))
        if (isTRUE(strict)) stop(msg, call. = FALSE)
        skipped <<- rbind(skipped, data.frame(method = method, reason = "execution_failed", stringsAsFactors = FALSE))
        messages <<- c(messages, msg)
        NULL
      }
    )

    if (is.null(mat)) next

    mat <- align_deconv_to_seurat(mat, seu = seu, method_name = method)
    if (isTRUE(normalize)) {
      mat <- normalize_deconv_rows(mat)
    }
    deconv_out[[method]] <- attach_method_metadata(
      mat,
      method_name = method,
      imported_from = "runner",
      normalized = isTRUE(normalize),
      original_file = NA_character_
    )
  }

  if (length(deconv_out) == 0L) {
    stop(
      "No methods produced runnable deconvolution outputs. Review `methods`, dependencies, and strict/use_python settings.",
      call. = FALSE
    )
  }

  list(
    seu = seu,
    deconv = deconv_out,
    methods_run = names(deconv_out),
    methods_skipped = skipped,
    messages = messages,
    sample_id = sample_id
  )
}

#' Run full AEGIS pipeline from raw inputs in one call
#'
#' Convenience wrapper that runs selected deconvolution methods first, then
#' passes standardized outputs into `run_aegis()`.
#'
#' @param seu A Seurat spatial object.
#' @param reference Reference input for runnable methods.
#' @param methods Character vector of method names.
#' @param markers Optional marker list.
#' @param strict Logical; strict behavior in `run_deconvolution()`.
#' @param use_python Logical; allow Python-backed methods.
#' @param normalize Logical; normalize method outputs.
#' @param ... Additional arguments forwarded to `run_deconvolution()`.
#'
#' @return An `aegis` object with downstream analyses completed.
#' @export
run_aegis_full <- function(
    seu,
    reference,
    methods = c("SPOTlight", "RCTD", "CARD"),
    markers = NULL,
    strict = TRUE,
    use_python = TRUE,
    normalize = TRUE,
    ...) {
  res <- run_deconvolution(
    seu = seu,
    reference = reference,
    methods = methods,
    normalize = normalize,
    strict = strict,
    use_python = use_python,
    ...
  )

  obj <- run_aegis(
    x = res$seu,
    deconv = res$deconv,
    markers = markers,
    do_marker = TRUE,
    do_spatial = TRUE,
    do_compare = TRUE,
    do_consensus = TRUE
  )

  obj$meta$deconvolution_run <- list(
    methods_run = res$methods_run,
    methods_skipped = res$methods_skipped,
    messages = res$messages
  )
  obj
}

#' Run SPOTlight deconvolution from R
#'
#' @param seu A Seurat spatial object.
#' @param reference Reference input accepted by backend workflow.
#' @param normalize Logical; row-normalize output.
#' @param strict Logical; strict parsing behavior.
#' @param ... Additional backend arguments.
#'
#' @return Spot-by-celltype numeric matrix.
#' @export
run_spotlight <- function(
    seu,
    reference,
    normalize = TRUE,
    strict = TRUE,
    ...) {
  run_r_native_method(
    seu = seu,
    reference = reference,
    method = "SPOTlight",
    package = "SPOTlight",
    function_candidates = c("SPOTlight", "runDeconvolution"),
    normalize = normalize,
    strict = strict,
    ...
  )
}

#' Run CARD deconvolution from R
#'
#' @param seu A Seurat spatial object.
#' @param reference Reference input accepted by backend workflow.
#' @param normalize Logical; row-normalize output.
#' @param strict Logical; strict parsing behavior.
#' @param ... Additional backend arguments.
#'
#' @return Spot-by-celltype numeric matrix.
#' @export
run_card <- function(
    seu,
    reference,
    normalize = TRUE,
    strict = TRUE,
    ...) {
  run_r_native_method(
    seu = seu,
    reference = reference,
    method = "CARD",
    package = "CARD",
    function_candidates = c("CARD_deconvolution", "CARD_deconv"),
    normalize = normalize,
    strict = strict,
    ...
  )
}

#' Run RCTD deconvolution from R
#'
#' @param seu A Seurat spatial object.
#' @param reference Reference input accepted by backend workflow.
#' @param normalize Logical; row-normalize output.
#' @param strict Logical; strict parsing behavior.
#' @param ... Additional backend arguments.
#'
#' @return Spot-by-celltype numeric matrix.
#' @export
run_rctd <- function(
    seu,
    reference,
    normalize = TRUE,
    strict = TRUE,
    ...) {
  run_r_native_method(
    seu = seu,
    reference = reference,
    method = "RCTD",
    package = "spacexr",
    function_candidates = c("run.RCTD", "run_RCTD", "runRCTD"),
    normalize = normalize,
    strict = strict,
    ...
  )
}

#' Run cell2location via reticulate (optional)
#'
#' @param seu A Seurat spatial object.
#' @param reference Reference input.
#' @param normalize Logical; row-normalize output.
#' @param strict Logical; strict parsing behavior.
#' @param ... Additional backend arguments.
#'
#' @return Spot-by-celltype numeric matrix.
#' @export
run_cell2location <- function(
    seu,
    reference,
    normalize = TRUE,
    strict = TRUE,
    ...) {
  run_python_method(
    seu = seu,
    reference = reference,
    method = "cell2location",
    modules = c("cell2location", "scvi"),
    normalize = normalize,
    strict = strict,
    ...
  )
}

#' Run DestVI via reticulate (optional)
#'
#' @param seu A Seurat spatial object.
#' @param reference Reference input.
#' @param normalize Logical; row-normalize output.
#' @param strict Logical; strict parsing behavior.
#' @param ... Additional backend arguments.
#'
#' @return Spot-by-celltype numeric matrix.
#' @export
run_destvi <- function(
    seu,
    reference,
    normalize = TRUE,
    strict = TRUE,
    ...) {
  run_python_method(
    seu = seu,
    reference = reference,
    method = "DestVI",
    modules = c("scvi"),
    normalize = normalize,
    strict = strict,
    ...
  )
}

#' Run Tangram via reticulate (optional)
#'
#' @param seu A Seurat spatial object.
#' @param reference Reference input.
#' @param normalize Logical; row-normalize output.
#' @param strict Logical; strict parsing behavior.
#' @param ... Additional backend arguments.
#'
#' @return Spot-by-celltype numeric matrix.
#' @export
run_tangram <- function(
    seu,
    reference,
    normalize = TRUE,
    strict = TRUE,
    ...) {
  run_python_method(
    seu = seu,
    reference = reference,
    method = "Tangram",
    modules = c("tangram"),
    normalize = normalize,
    strict = strict,
    ...
  )
}

#' Run stereoscope via reticulate (optional)
#'
#' @param seu A Seurat spatial object.
#' @param reference Reference input.
#' @param normalize Logical; row-normalize output.
#' @param strict Logical; strict parsing behavior.
#' @param ... Additional backend arguments.
#'
#' @return Spot-by-celltype numeric matrix.
#' @export
run_stereoscope <- function(
    seu,
    reference,
    normalize = TRUE,
    strict = TRUE,
    ...) {
  run_python_method(
    seu = seu,
    reference = reference,
    method = "stereoscope",
    modules = c("scvi"),
    normalize = normalize,
    strict = strict,
    ...
  )
}

#' @keywords internal
resolve_method_names <- function(methods, registry = get_supported_methods(), unique_out = TRUE) {
  normalize_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))

  aliases <- setNames(registry$method_name, normalize_key(registry$method_name))
  aliases[[normalize_key("cell2loc")]] <- "cell2location"
  aliases[[normalize_key("spotlight")]] <- "SPOTlight"
  aliases[[normalize_key("stdeconv")]] <- "STdeconvolve"

  resolved <- vapply(methods, function(m) {
    key <- normalize_key(m)
    if (!key %in% names(aliases)) return(NA_character_)
    aliases[[key]]
  }, character(1))

  bad <- methods[is.na(resolved)]
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "Unsupported method(s): %s. Use get_supported_methods() to see valid names.",
        paste(unique(bad), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  out <- unname(resolved)
  if (isTRUE(unique_out)) out <- unique(out)
  out
}

#' @keywords internal
pick_method_runner <- function(method, registry, runner_overrides = NULL) {
  if (!is.null(runner_overrides) && method %in% names(runner_overrides)) {
    fn <- runner_overrides[[method]]
    if (!is.function(fn)) {
      stop(sprintf("runner_overrides[['%s']] must be a function.", method), call. = FALSE)
    }
    return(fn)
  }

  reg_row <- registry[registry$method_name == method, , drop = FALSE]
  runner_name <- reg_row$runner_function[[1L]]
  if (is.na(runner_name) || !nzchar(runner_name)) {
    stop(sprintf("Method '%s' has no runner function registered.", method), call. = FALSE)
  }

  fn <- get0(runner_name, envir = asNamespace("AEGIS"), inherits = FALSE)
  if (!is.function(fn)) {
    stop(sprintf("Registered runner '%s' for method '%s' was not found.", runner_name, method), call. = FALSE)
  }
  fn
}

#' @keywords internal
validate_reference_input <- function(reference, method = "method") {
  if (missing(reference) || is.null(reference)) {
    stop(sprintf("Method '%s' requires a non-NULL `reference` input.", method), call. = FALSE)
  }
  ok <- inherits(reference, "Seurat") ||
    is.matrix(reference) ||
    is.data.frame(reference) ||
    (is.list(reference) && length(reference) > 0L)
  if (!ok) {
    stop(
      sprintf(
        "Method '%s': unsupported `reference` type '%s'. Expected Seurat, matrix/data.frame, or non-empty list.",
        method,
        paste(class(reference), collapse = ",")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
run_r_native_method <- function(
    seu,
    reference,
    method,
    package,
    function_candidates,
    normalize = TRUE,
    strict = TRUE,
    ...) {
  assert_is_seurat(seu, "seu")
  validate_reference_input(reference, method = method)

  dots <- list(...)
  if ("result" %in% names(dots) && !is.null(dots$result)) {
    return(standardize_runner_output(
      result = dots$result,
      method = method,
      seurat_spots = colnames(seu),
      normalize = normalize,
      strict = strict,
      spot_col = dots$spot_col %||% NULL
    ))
  }

  if ("run_fn" %in% names(dots) && is.function(dots$run_fn)) {
    dots_for_fn <- dots
    dots_for_fn[c("run_fn", "result", "backend_args", "spot_col")] <- NULL
    raw <- do.call(dots$run_fn, c(list(seu = seu, reference = reference), dots_for_fn))
    return(standardize_runner_output(
      result = raw,
      method = method,
      seurat_spots = colnames(seu),
      normalize = normalize,
      strict = strict,
      spot_col = dots$spot_col %||% NULL
    ))
  }

  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      sprintf(
        "Method '%s' requires R package '%s'. Install it to run directly, or import exported outputs with %s().",
        method,
        package,
        adapter_for_method(method)
      ),
      call. = FALSE
    )
  }

  backend_fn <- NULL
  ns <- asNamespace(package)
  for (fn_name in function_candidates) {
    fn <- get0(fn_name, envir = ns, inherits = FALSE)
    if (is.function(fn)) {
      backend_fn <- fn
      break
    }
  }

  if (is.null(backend_fn)) {
    stop(
      sprintf(
        "Method '%s': could not find an expected backend function in package '%s' (%s). Provide `run_fn` explicitly or import exported results via %s().",
        method,
        package,
        paste(function_candidates, collapse = ", "),
        adapter_for_method(method)
      ),
      call. = FALSE
    )
  }

  backend_args <- dots$backend_args
  if (is.null(backend_args)) backend_args <- list()
  if (!is.list(backend_args)) {
    stop("`backend_args` must be NULL or a list.", call. = FALSE)
  }

  fm <- names(formals(backend_fn))
  if (!is.null(fm)) {
    if ("seu" %in% fm && !("seu" %in% names(backend_args))) backend_args$seu <- seu
    if ("reference" %in% fm && !("reference" %in% names(backend_args))) backend_args$reference <- reference
    if ("spatial" %in% fm && !("spatial" %in% names(backend_args))) backend_args$spatial <- seu
    if ("query" %in% fm && !("query" %in% names(backend_args))) backend_args$query <- seu
    if ("x" %in% fm && !("x" %in% names(backend_args))) backend_args$x <- reference
    if ("y" %in% fm && !("y" %in% names(backend_args))) backend_args$y <- seu
  }

  raw <- tryCatch(
    do.call(backend_fn, backend_args),
    error = function(e) {
      stop(
        sprintf(
          "Method '%s' backend call failed (%s). Provide explicit `backend_args`/`run_fn`, or import exported outputs via %s().",
          method,
          conditionMessage(e),
          adapter_for_method(method)
        ),
        call. = FALSE
      )
    }
  )

  standardize_runner_output(
    result = raw,
    method = method,
    seurat_spots = colnames(seu),
    normalize = normalize,
    strict = strict,
    spot_col = dots$spot_col %||% NULL
  )
}

#' @keywords internal
run_python_method <- function(
    seu,
    reference,
    method,
    modules,
    normalize = TRUE,
    strict = TRUE,
    ...) {
  assert_is_seurat(seu, "seu")
  validate_reference_input(reference, method = method)

  dots <- list(...)
  if ("result" %in% names(dots) && !is.null(dots$result)) {
    return(standardize_runner_output(
      result = dots$result,
      method = method,
      seurat_spots = colnames(seu),
      normalize = normalize,
      strict = strict,
      spot_col = dots$spot_col %||% NULL
    ))
  }

  ready <- check_python_method_ready(method = method, modules = modules)
  if (!isTRUE(ready$ready)) {
    stop(ready$message, call. = FALSE)
  }

  if (!("run_fn" %in% names(dots)) || !is.function(dots$run_fn)) {
    stop(
      sprintf(
        "Method '%s' Python wrapper is available, but an executable default pipeline is not guaranteed in this environment. Provide `run_fn` (reticulate-backed) or import exported outputs via %s().",
        method,
        adapter_for_method(method)
      ),
      call. = FALSE
    )
  }

  dots_for_fn <- dots
  dots_for_fn[c("run_fn", "result", "spot_col")] <- NULL
  raw <- do.call(dots$run_fn, c(list(seu = seu, reference = reference), dots_for_fn))
  standardize_runner_output(
    result = raw,
    method = method,
    seurat_spots = colnames(seu),
    normalize = normalize,
    strict = strict,
    spot_col = dots$spot_col %||% NULL
  )
}

#' @keywords internal
check_python_method_ready <- function(method, modules = NULL) {
  if (is.null(modules)) {
    reg <- get_supported_methods()
    row <- reg[reg$method_name == method, , drop = FALSE]
    dep <- row$backend_dependency[[1L]] %||% ""
    modules <- trimws(unlist(strsplit(dep, ",", fixed = TRUE)))
    modules <- modules[nzchar(modules)]
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(list(
      method = method,
      ready = FALSE,
      missing_modules = modules,
      message = sprintf(
        "Method '%s' requires Python via `reticulate`, but package 'reticulate' is not installed. Install reticulate or import exported outputs via %s().",
        method,
        adapter_for_method(method)
      )
    ))
  }

  if (length(modules) == 0L) {
    return(list(method = method, ready = TRUE, missing_modules = character(0), message = "OK"))
  }

  missing <- modules[!vapply(modules, reticulate::py_module_available, logical(1))]
  if (length(missing) > 0L) {
    return(list(
      method = method,
      ready = FALSE,
      missing_modules = missing,
      message = sprintf(
        "Method '%s' is not Python-ready. Missing module(s): %s. Configure Python/reticulate env, or import exported outputs via %s().",
        method,
        paste(missing, collapse = ", "),
        adapter_for_method(method)
      )
    ))
  }

  list(method = method, ready = TRUE, missing_modules = character(0), message = "OK")
}

#' @keywords internal
standardize_runner_output <- function(
    result,
    method,
    seurat_spots,
    normalize = TRUE,
    strict = TRUE,
    spot_col = NULL) {
  obj <- result
  if (!(is.matrix(obj) || is.data.frame(obj))) {
    obj <- extract_matrix_like_from_object(obj, type = "auto", method = method, strict = strict)
  }

  mat <- standardize_deconv_matrix(obj = obj, spot_col = spot_col, strict = strict, method = method)
  mat <- maybe_transpose_deconv_matrix(mat, strict = strict, method = method, seurat_spots = seurat_spots)
  mat <- align_matrix_to_spots(mat, spots = seurat_spots, method_name = method)
  if (isTRUE(normalize)) {
    mat <- normalize_deconv_rows(mat)
  }

  assert_numeric_deconv(mat, method = method)
}

#' @keywords internal
adapter_for_method <- function(method) {
  reg <- get_supported_methods()
  row <- reg[reg$method_name == method, , drop = FALSE]
  if (nrow(row) == 0L) return("read_deconv_table")
  row$adapter_reader[[1L]] %||% "read_deconv_table"
}
