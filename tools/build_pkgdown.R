#!/usr/bin/env Rscript

# Local pkgdown helper. If Pandoc is unavailable, build non-article pages.

if (!requireNamespace("pkgdown", quietly = TRUE)) {
  stop("pkgdown is not installed. Install with install.packages('pkgdown').", call. = FALSE)
}

if (rmarkdown::pandoc_available()) {
  pkgdown::build_site(new_process = FALSE, install = TRUE)
} else {
  message("Pandoc not found; skip local pkgdown build. Use GitHub Actions workflow 'pkgdown' for full site generation.")
}
