# Changelog

## AEGIS 0.1.0

- Initial MVP release.
- Added Seurat-first data loading for Human Lymph Node 10x Visium
  assets.
- Added simulated deconvolution, audit, comparison, plotting, and report
  pipeline.
- Added example data, vignettes, and test coverage for core workflow.
- Updated README with enlarged logo, devtools installation, citation
  text, and repository link.
- Added tutorial expansion, spatial slice example figure, and
  pkgdown/R-CMD-check workflows.
- Added unified
  [`run_aegis()`](https://jameswu7.github.io/AEGIS/reference/run_aegis.md)
  one-call workflow for single-sample and multi-sample usage.
- Updated tutorials to include multi-sample end-to-end examples.
- P7 pre-release hardening: added integration workflow tests
  (simulated/import/multi-sample), docs consistency checks, and CI smoke
  workflow.
- Added robust tutorial link fallbacks and a reproducible README figure
  regeneration script.
- Improved importer error messages for ambiguous spot IDs and missing
  numeric cell-type columns.
