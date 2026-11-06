# altiro 0.6.3

* Documentation and messages call the classification "classification rules" instead of
  "decision tree", to avoid confusion with the machine-learning algorithm. The description of
  `unclassified` features changes accordingly; labels and results are unchanged.
* Manual: new subsection on the significance threshold `alpha` (Python, R and command line,
  re-classification without refitting, effect on primed calls).

# altiro 0.6.2

* `altiro classify` (command line) gains `--lfc-col`, `--padj-col` and `--se-col`, so contrast tables
  from edgeR (`logFC`, `FDR`) or limma (`logFC`, `adj.P.Val`) are read directly.
* Manual: the section on precomputed contrasts now shows complete edgeR and limma-voom workflows in R,
  the Python and command-line equivalents, and the column names of each tool.

# altiro 0.6.1

* Manual: new subsection on the optional length correction of `merge_ortholog_counts()`
  (how to pass `lengths_a`/`lengths_b`, integer output, caveats for DESeq2).
* Tests now check that length-normalised matrices stay integer in Python and R.

# altiro 0.6.0

* The package is renamed from `orthoresponse` to `altiro` (ALTIRO, Analysis of Latent Transcriptional Induction and Response Outcomes).
  R: `orthoresponse()` is now `altiro()` and `save_orthoresponse()` is `save_altiro()`.
  Python: `import altiro`, `run_orthoresponse()` is now `run_altiro()`, `OrthoResponseResult` is `AltiroResult`,
  and the command line tool is `altiro`. Results, options and file formats are unchanged.
* The GitHub repository moves to `Cold7/altiro`.

# altiro 0.5.1

* `altiro()` / `run_altiro()` gain `group_names` to set display names of groups A and B independently of the sample-table values; the CLI option `run --group-names` now works.
* Manual: new section on customizing names and labels.

# altiro 0.5.0

* The stress-ready model is now called **primed** (`primed_1` ... `primed_4`), a neutral name for any treatment, organism or group. The classification is unchanged.
* `classify_responses(label_style = "stress_ready")` (and the Python/CLI equivalents) returns the old names.
* `primed_asymmetry()` columns are now `primed` and `fraction_primed`; `plot_model_patterns()` uses the row label "Primed" and the key `primed` in `model_labels` (the old key `stress_ready` is still accepted).

# altiro 0.4.1

* Python: `print(res)` shows the same summary as `print(res)` in R, including the stress-ready asymmetry line with its Fisher p-value; the command line prints it too. New `format_asymmetry()`.

# altiro 0.4.0

* `save_altiro()` (R) and `AltiroResult.save()` (Python) also save the figures `model_patterns` and `model_summary` (PNG and PDF) by default; options `plots`, `condition_labels`, `model_labels`, `group_names`, `title`, `formats`.
* New `save_plots()` in both languages; Python CLI options `--condition-labels`, `--plot-formats`, `--no-plots`.
* `plot_summary()` (R) leaves room for the value labels.

# altiro 0.3.0

* New bundled example `ovarian_brequinar` (GEO GSE283381): OVCAR5 cisplatin-sensitive vs -resistant cells, vehicle vs brequinar, raw HTSeq counts deposited by the authors. It replaces the GSE232890 example.
* README documents where the example files are and how to load them from the package or from file paths.

# altiro 0.2.0

* Example data bundled with the package: rice salt stress, GEO GSE167342 (`load_example()`, `example_files()`, `export_example()`, `list_examples()`; Python CLI: `altiro example`).

# altiro 0.1.0

* First release: decision tree of Contreras-Riquelme et al. (2025) (`preset = "paper"`), strict options (missing values, basal difference, equivalence test), DESeq2 wrapper, reciprocal best hits, ortholog count merging, summaries and plots.
* `plot_model_patterns()` accepts custom `condition_labels` (default `c("C", "S")`) and `title`.
* `plot_model_patterns()` accepts custom row labels (`model_labels`) and draws shared and opposite responses in the first and fourth columns.
