#!/usr/bin/env Rscript

# Local check helper. Use GitHub Actions for full check with Pandoc.

args <- c("--no-manual", "--ignore-vignettes")
res <- system2("R", c("CMD", "check", ".", args), stdout = TRUE, stderr = TRUE)
cat(paste(res, collapse = "\n"), "\n")
