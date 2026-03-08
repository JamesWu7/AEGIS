# Package index

## Data input

- [`load_10x_lymphnode()`](https://jameswu7.github.io/AEGIS/reference/load_10x_lymphnode.md)
  : Load the 10x Human Lymph Node Visium dataset as a Seurat object
- [`aegis_example`](https://jameswu7.github.io/AEGIS/reference/aegis_example.md)
  : Example Seurat object derived from Human Lymph Node data
- [`aegis_default_markers()`](https://jameswu7.github.io/AEGIS/reference/aegis_default_markers.md)
  : Default marker list used by AEGIS MVP

## Core object and simulation

- [`as_aegis()`](https://jameswu7.github.io/AEGIS/reference/as_aegis.md)
  : Construct an AEGIS object
- [`is_aegis()`](https://jameswu7.github.io/AEGIS/reference/is_aegis.md)
  : AEGIS analysis object
- [`simulate_deconv_results()`](https://jameswu7.github.io/AEGIS/reference/simulate_deconv_results.md)
  : Simulate deconvolution results for a Seurat spatial object

## Audit and comparison

- [`audit_basic()`](https://jameswu7.github.io/AEGIS/reference/audit_basic.md)
  : Run basic quality audits on deconvolution matrices
- [`audit_marker()`](https://jameswu7.github.io/AEGIS/reference/audit_marker.md)
  : Run marker-based concordance audit
- [`audit_spatial()`](https://jameswu7.github.io/AEGIS/reference/audit_spatial.md)
  : Run spatial smoothness audit
- [`compare_methods()`](https://jameswu7.github.io/AEGIS/reference/compare_methods.md)
  : Compare deconvolution methods
- [`compute_consensus()`](https://jameswu7.github.io/AEGIS/reference/compute_consensus.md)
  : Compute consensus deconvolution profile

## Visualization and reporting

- [`plot_audit()`](https://jameswu7.github.io/AEGIS/reference/plot_audit.md)
  : Plot audit outputs
- [`plot_compare()`](https://jameswu7.github.io/AEGIS/reference/plot_compare.md)
  : Plot method comparison and consensus outputs
- [`render_report()`](https://jameswu7.github.io/AEGIS/reference/render_report.md)
  : Render an AEGIS HTML report
- [`build_aegis_example_data()`](https://jameswu7.github.io/AEGIS/reference/build_aegis_example_data.md)
  : Build and save example data from raw Human Lymph Node files
