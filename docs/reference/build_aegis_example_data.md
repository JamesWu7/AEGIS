# Build and save example data from raw Human Lymph Node files

Build and save example data from raw Human Lymph Node files

## Usage

``` r
build_aegis_example_data(
  raw_dir = ".",
  output_file = file.path("data", "aegis_example.rda"),
  object_name = "aegis_example",
  max_spots = 1200L,
  seed = 42
)
```

## Arguments

- raw_dir:

  Directory containing authoritative raw files.

- output_file:

  Path to `.rda` output file.

- object_name:

  Name of the object saved in `.rda`.

- max_spots:

  Maximum number of spots to keep for package-size control.

- seed:

  Random seed for reproducible spot subsetting.

## Value

Invisibly returns the saved Seurat object.
