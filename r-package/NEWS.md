# orthoresponse 0.5.0

* The stress-ready model is now called **primed** (`primed_1` ... `primed_4`), a neutral name for any treatment, organism or group. The classification is unchanged.
* `classify_responses(label_style = "stress_ready")` (and the Python/CLI equivalents) returns the old names.
* `primed_asymmetry()` columns are now `primed` and `fraction_primed`; `plot_model_patterns()` uses the row label "Primed" and the key `primed` in `model_labels` (the old key `stress_ready` is still accepted).

# orthoresponse 0.4.1

* Python: `print(res)` shows the same summary as `print(res)` in R, including the stress-ready asymmetry line with its Fisher p-value; the command line prints it too. New `format_asymmetry()`.

# orthoresponse 0.4.0

* `save_orthoresponse()` (R) and `OrthoResponseResult.save()` (Python) also save the figures `model_patterns` and `model_summary` (PNG and PDF) by default; options `plots`, `condition_labels`, `model_labels`, `group_names`, `title`, `formats`.
* New `save_plots()` in both languages; Python CLI options `--condition-labels`, `--plot-formats`, `--no-plots`.
* `plot_summary()` (R) leaves room for the value labels.

# orthoresponse 0.3.0

* New bundled example `ovarian_brequinar` (GEO GSE283381): OVCAR5 cisplatin-sensitive vs -resistant cells, vehicle vs brequinar, raw HTSeq counts deposited by the authors. It replaces the GSE232890 example.
* README documents where the example files are and how to load them from the package or from file paths.

# orthoresponse 0.2.0

* Example data bundled with the package: rice salt stress, GEO GSE167342 (`load_example()`, `example_files()`, `export_example()`, `list_examples()`; Python CLI: `orthoresponse example`).

# orthoresponse 0.1.0

* First release: decision tree of Contreras-Riquelme et al. (2025) (`preset = "paper"`), strict options (missing values, basal difference, equivalence test), DESeq2 wrapper, reciprocal best hits, ortholog count merging, summaries and plots.
* `plot_model_patterns()` accepts custom `condition_labels` (default `c("C", "S")`) and `title`.
* `plot_model_patterns()` accepts custom row labels (`model_labels`) and draws shared and opposite responses in the first and fourth columns.
