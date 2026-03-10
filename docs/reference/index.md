# Package index

## Data input

- [`load_10x_lymphnode()`](https://jameswu7.github.io/AEGIS/reference/load_10x_lymphnode.md)
  : Load the 10x Human Lymph Node Visium dataset as a Seurat object
- [`load_10x_spatial_set()`](https://jameswu7.github.io/AEGIS/reference/load_10x_spatial_set.md)
  : Load a set of 10x spatial datasets
- [`merge_spatial_seurat_list()`](https://jameswu7.github.io/AEGIS/reference/merge_spatial_seurat_list.md)
  : Merge a list of spatial Seurat objects
- [`aegis_example`](https://jameswu7.github.io/AEGIS/reference/aegis_example.md)
  : Example Seurat object derived from Human Lymph Node data
- [`aegis_default_markers()`](https://jameswu7.github.io/AEGIS/reference/aegis_default_markers.md)
  : Default marker list used by AEGIS MVP
- [`align_deconv_to_seurat()`](https://jameswu7.github.io/AEGIS/reference/align_deconv_to_seurat.md)
  : Align a deconvolution matrix to Seurat spot order
- [`read_rctd()`](https://jameswu7.github.io/AEGIS/reference/read_rctd.md)
  : Read RCTD deconvolution outputs
- [`read_spotlight()`](https://jameswu7.github.io/AEGIS/reference/read_spotlight.md)
  : Read SPOTlight deconvolution outputs
- [`read_cell2location()`](https://jameswu7.github.io/AEGIS/reference/read_cell2location.md)
  : Read cell2location deconvolution outputs

## Core object and simulation

- [`run_aegis()`](https://jameswu7.github.io/AEGIS/reference/run_aegis.md)
  : Run the standard AEGIS workflow with one function
- [`as_aegis()`](https://jameswu7.github.io/AEGIS/reference/as_aegis.md)
  : Construct an AEGIS object
- [`as_aegis_multi()`](https://jameswu7.github.io/AEGIS/reference/as_aegis_multi.md)
  : Construct a multi-sample AEGIS object
- [`is_aegis()`](https://jameswu7.github.io/AEGIS/reference/is_aegis.md)
  : AEGIS analysis object
- [`is_aegis_multi()`](https://jameswu7.github.io/AEGIS/reference/is_aegis_multi.md)
  : Check if object is an AEGIS multi-sample object
- [`split_aegis_by_sample()`](https://jameswu7.github.io/AEGIS/reference/split_aegis_by_sample.md)
  : Split an AEGIS object by sample
- [`summarize_by_sample()`](https://jameswu7.github.io/AEGIS/reference/summarize_by_sample.md)
  : Summarize AEGIS results by sample
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
- [`render_report_batch()`](https://jameswu7.github.io/AEGIS/reference/render_report_batch.md)
  : Render one report per sample
- [`build_aegis_example_data()`](https://jameswu7.github.io/AEGIS/reference/build_aegis_example_data.md)
  : Build and save example data from raw Human Lymph Node files
