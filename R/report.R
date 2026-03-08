#' Render an AEGIS HTML report
#'
#' Renders a standalone HTML report from the package template.
#'
#' @param x An `aegis` object.
#' @param output_file Output HTML path.
#'
#' @return Path to rendered HTML report.
#' @export
render_report <- function(x, output_file = "aegis_report.html") {
  assert_is_aegis(x)
  template <- resolve_template()
  output_file <- normalizePath(output_file, winslash = "/", mustWork = FALSE)
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  render_env <- new.env(parent = baseenv())
  render_env$params <- list(aegis_obj = x)
  render_env$audit_basic <- audit_basic
  render_env$audit_marker <- audit_marker
  render_env$audit_spatial <- audit_spatial
  render_env$compare_methods <- compare_methods
  render_env$compute_consensus <- compute_consensus
  render_env$plot_audit <- plot_audit
  render_env$plot_compare <- plot_compare
  render_env$head <- utils::head
  render_env$sessionInfo <- utils::sessionInfo

  if (rmarkdown::pandoc_available()) {
    rendered <- rmarkdown::render(
      input = template,
      output_file = basename(output_file),
      output_dir = out_dir,
      params = list(aegis_obj = x),
      envir = render_env,
      quiet = TRUE
    )
    return(normalizePath(rendered, winslash = "/", mustWork = TRUE))
  }

  if (!requireNamespace("markdown", quietly = TRUE)) {
    stop(
      "Pandoc is unavailable and package 'markdown' is not installed. Install Pandoc or install 'markdown' for fallback rendering.",
      call. = FALSE
    )
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)

  md_file <- tempfile(pattern = "aegis_report_", fileext = ".md", tmpdir = out_dir)
  knitr::knit(
    input = template,
    output = md_file,
    envir = render_env,
    quiet = TRUE
  )
  markdown::markdownToHTML(
    file = md_file,
    output = basename(output_file),
    options = c("use_xhtml", "base64_images")
  )

  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}
