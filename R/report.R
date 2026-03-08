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

  if (rmarkdown::pandoc_available()) {
    rendered <- rmarkdown::render(
      input = template,
      output_file = basename(output_file),
      output_dir = out_dir,
      params = list(aegis_obj = x),
      envir = new.env(parent = baseenv()),
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

  env <- new.env(parent = baseenv())
  env$params <- list(aegis_obj = x)
  md_file <- tempfile(fileext = ".md")
  knitr::knit(input = template, output = md_file, envir = env, quiet = TRUE)
  markdown::markdownToHTML(
    file = md_file,
    output = output_file,
    options = c("use_xhtml", "base64_images")
  )

  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}
