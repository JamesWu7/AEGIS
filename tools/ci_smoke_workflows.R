#!/usr/bin/env Rscript

required_pkgs <- c("devtools")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop(sprintf("Missing required package(s): %s", paste(missing, collapse = ", ")), call. = FALSE)
}

devtools::load_all(".", quiet = TRUE)

data("aegis_example", package = "AEGIS")
markers <- readRDS(system.file("extdata", "marker_list.rds", package = "AEGIS"))

# Single-sample simulated smoke.
deconv <- simulate_deconv_results(aegis_example, seed = 777)
obj <- run_aegis(aegis_example, deconv = deconv, markers = markers)
stopifnot(inherits(obj, "aegis"))
stopifnot(inherits(plot_audit(obj, "dominance"), "ggplot"))
stopifnot(inherits(plot_compare(obj, "heatmap"), "ggplot"))
out <- tempfile(fileext = ".html")
report <- render_report(obj, output_file = out)
stopifnot(file.exists(report))

# Import workflow smoke with fixture-like tables.
spots <- colnames(aegis_example)[1:6]
seu_small <- suppressWarnings(aegis_example[, spots])

tmp_rctd <- tempfile(fileext = ".csv")
utils::write.csv(
  data.frame(
    barcode = spots,
    B_cell = c(0.4, 0.2, 0.1, 0.3, 0.2, 0.5),
    T_cell = c(0.4, 0.5, 0.6, 0.4, 0.6, 0.3),
    Myeloid = c(0.2, 0.3, 0.3, 0.3, 0.2, 0.2),
    check.names = FALSE
  ),
  tmp_rctd,
  row.names = FALSE
)
rctd <- read_rctd(tmp_rctd)
obj2 <- as_aegis(seu_small, deconv = list(RCTD = rctd), markers = markers)
obj2 <- audit_basic(obj2)
stopifnot(inherits(obj2, "aegis"))
stopifnot(!is.null(obj2$audit$basic$summary))

message("AEGIS CI smoke workflows succeeded.")
