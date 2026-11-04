# orthoresponse 0.2.0

* Example data bundled with the package: rice salt stress, GEO GSE167342 (`load_example()`, `example_files()`, `export_example()`, `list_examples()`; Python CLI: `orthoresponse example`).

# orthoresponse 0.1.0

* First release: decision tree of Contreras-Riquelme et al. (2025) (`preset = "paper"`), strict options (missing values, basal difference, equivalence test), DESeq2 wrapper, reciprocal best hits, ortholog count merging, summaries and plots.
* `plot_model_patterns()` accepts custom `condition_labels` (default `c("C", "S")`) and `title`.
* `plot_model_patterns()` accepts custom row labels (`model_labels`) and draws shared and opposite responses in the first and fourth columns.
