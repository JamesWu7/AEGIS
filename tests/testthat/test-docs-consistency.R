find_repo_root <- function() {
  inst_pkg <- system.file(package = "AEGIS")
  candidates <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", ".."),
    if (nzchar(inst_pkg)) inst_pkg else NA_character_,
    if (nzchar(inst_pkg)) file.path(inst_pkg, "..", "..") else NA_character_
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  for (cand in candidates) {
    if (!dir.exists(cand)) next
    cand_norm <- normalizePath(cand, winslash = "/", mustWork = TRUE)
    has_desc <- file.exists(file.path(cand_norm, "DESCRIPTION"))
    has_markers <- (
      file.exists(file.path(cand_norm, "README.md")) ||
      dir.exists(file.path(cand_norm, "vignettes")) ||
      dir.exists(file.path(cand_norm, "docs"))
    )
    if (has_desc && has_markers) {
      return(cand_norm)
    }
  }
  NA_character_
}

test_that("README tutorial references and figure assets exist", {
  repo_root <- find_repo_root()
  if (is.na(repo_root)) {
    testthat::skip("Repository root is not available in this test context.")
  }
  readme_path <- file.path(repo_root, "README.md")

  if (!file.exists(readme_path)) {
    # Installed-package context (R CMD check): README is typically unavailable.
    testthat::skip("README.md is not available in installed-package test context.")
  }

  readme <- tryCatch(
    paste(readLines(readme_path, warn = FALSE), collapse = "\n"),
    error = function(e) NA_character_
  )
  if (is.na(readme)) {
    testthat::skip("README.md could not be read in this test context.")
  }

  expect_match(readme, "Tutorials")
  expect_match(readme, "AEGIS-overview")
  expect_match(readme, "AEGIS-complete-tutorial")
  expect_match(readme, "\\[source\\]\\(vignettes/AEGIS-overview.Rmd\\)")
  expect_match(readme, "https://jameswu7.github.io/AEGIS/articles/AEGIS-overview.html")
  expect_match(readme, "https://htmlpreview.github.io/\\?https://github.com/JamesWu7/AEGIS/blob/main/docs/articles/AEGIS-overview.html")

  # Ensure README does not regress to bare vignette path inventory.
  expect_false(grepl("^\\s*-\\s*vignettes/", readme, perl = TRUE))

  fig_paths <- c(
    file.path(repo_root, "inst", "assets", "figures", "readme-slice.png"),
    file.path(repo_root, "inst", "assets", "figures", "readme-dominance.png")
  )
  for (p in fig_paths) {
    if (file.exists(p)) {
      expect_true(file.exists(p), info = sprintf("Missing README asset: %s", p))
      expect_gt(file.info(p)$size, 0)
    }
  }

  regen_script <- file.path(repo_root, "tools", "regenerate_readme_figures.R")
  if (file.exists(regen_script)) {
    expect_true(file.exists(regen_script))
  }
})

test_that("README-listed key API functions are exported", {
  ns <- getNamespaceExports("AEGIS")
  must_have <- c(
    "load_10x_lymphnode",
    "load_10x_spatial_set",
    "simulate_deconv_results",
    "read_deconv_table",
    "read_rctd",
    "read_spotlight",
    "read_cell2location",
    "read_card",
    "read_spatialdwls",
    "read_stereoscope",
    "read_destvi",
    "read_tangram",
    "read_stdeconvolve",
    "read_dstg",
    "read_stride",
    "as_aegis",
    "run_aegis",
    "audit_basic",
    "audit_marker",
    "audit_spatial",
    "compare_methods",
    "compute_consensus",
    "plot_audit",
    "plot_compare",
    "render_report",
    "render_report_batch",
    "summarize_by_sample"
  )
  expect_true(all(must_have %in% ns))
})

test_that("vignette and docs tutorial files are present when available", {
  repo_root <- find_repo_root()
  if (is.na(repo_root)) {
    testthat::skip("Could not resolve repository root in this test context.")
  }

  vignettes <- c(
    "vignettes/AEGIS-overview.Rmd",
    "vignettes/AEGIS-complete-tutorial.Rmd"
  )
  docs_articles <- c(
    "docs/articles/AEGIS-overview.html",
    "docs/articles/AEGIS-complete-tutorial.html"
  )

  has_source_vignettes <- all(file.exists(file.path(repo_root, "vignettes", c(
    "AEGIS-overview.Rmd",
    "AEGIS-complete-tutorial.Rmd"
  ))))

  if (isTRUE(has_source_vignettes)) {
    for (p in vignettes) {
      abs <- file.path(repo_root, p)
      expect_true(file.exists(abs), info = sprintf("Missing tutorial artifact: %s", p))
    }
  }

  docs_dir <- file.path(repo_root, "docs", "articles")
  if (dir.exists(docs_dir)) {
    for (p in docs_articles) {
      abs <- file.path(repo_root, p)
      expect_true(file.exists(abs), info = sprintf("Missing tutorial artifact: %s", p))
    }
  } else {
    # Installed package context: look for built vignette artifacts.
    inst_doc <- system.file("doc", package = "AEGIS")
    if (nzchar(inst_doc) && dir.exists(inst_doc)) {
      for (nm in c("AEGIS-overview.html", "AEGIS-complete-tutorial.html")) {
        expect_true(file.exists(file.path(inst_doc, nm)), info = sprintf("Missing installed tutorial artifact: %s", nm))
      }
    } else {
      testthat::skip("docs/articles or installed doc artifacts are not available in this context.")
    }
  }
})
