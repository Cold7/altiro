# orthoresponse

Classify genes or ortholog pairs into **expression response models** in a
2 x 2 design: two groups (species, genotypes, cultivars, cell lines) measured
under a control and a stress condition.

![Expression models](docs/figures/model_patterns.png)

| model | meaning |
|---|---|
| **stress-ready** (1-4) | one group responds and reaches the *constitutive* level of the other group, which is therefore primed |
| **shared** (1-2) | both groups respond in the same direction |
| **unique** (1-4) | only one group responds |
| **opposite** (1-2) | both respond, in opposite directions |
| **no response** (1-3) | neither responds (basal level higher, lower or equal) |

The default rules reproduce exactly the decision tree of
Contreras-Riquelme *et al.* (2025) *PLoS One* 20(5): e0324724
([doi:10.1371/journal.pone.0324724](https://doi.org/10.1371/journal.pone.0324724)).
Optional stricter rules handle missing statistics, require a basal difference
for stress-ready calls, and can replace "not significantly different" by a
formal equivalence test.

The same functionality is available in **Python** and **R** with the same
function names; both give identical classifications for identical inputs.

## Installation

**Python** (>= 3.9)

```bash
pip install "orthoresponse[all] @ git+https://github.com/Cold7/orthoresponse.git#subdirectory=python"
```

**R** (>= 4.1)

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("DESeq2")
remotes::install_github("Cold7/orthoresponse", subdir = "r-package", build_vignettes = TRUE)
```

## Example data

Two real data sets are included, both with raw counts deposited by the authors
in GEO and 3 replicates per group:

| name | data set | group A vs group B | control vs stress | size |
|---|---|---|---|---|
| `rice_salt` | GEO GSE167342, rice leaves | PL12 (salt-sensitive) vs IL22 (*Saltol* introgression line, tolerant) | Control vs Salt (24 h, 80 mM NaCl) | 36,844 genes x 12 samples |
| `ovarian_brequinar` | GEO GSE283381, OVCAR5 ovarian cancer cells | CisS (cisplatin-sensitive) vs CisR (cisplatin-resistant) | Vehicle vs Brequinar (DHODH inhibitor) | 59,453 genes x 12 samples |

### Where the files are

They are bundled inside the packages, not in `examples/`:

| file | in this repository (Python copy; the R copy is identical) | after installing |
|---|---|---|
| rice counts | `python/src/orthoresponse/data/rice_GSE167342_counts.tsv.gz` | Python: `orr.example_files("rice_salt")["counts"]`; R: `example_files("rice_salt")[["counts"]]` |
| rice samples | `python/src/orthoresponse/data/rice_GSE167342_samples.tsv` | Python: `orr.example_files("rice_salt")["samples"]`; R: `example_files("rice_salt")[["samples"]]` |
| ovarian counts | `python/src/orthoresponse/data/ovarian_GSE283381_counts.tsv.gz` | Python: `orr.example_files("ovarian_brequinar")["counts"]`; R: `example_files("ovarian_brequinar")[["counts"]]` |
| ovarian samples | `python/src/orthoresponse/data/ovarian_GSE283381_samples.tsv` | Python: `orr.example_files("ovarian_brequinar")["samples"]`; R: `example_files("ovarian_brequinar")[["samples"]]` |

The R copies are in `r-package/inst/extdata/` with the same file names. Count
files are tab-separated and gzip-compressed (genes in rows, samples in
columns); pandas and R read them directly. Each sample table maps every column
to its group and condition, for example (rice):

```text
sample  line  treatment  geo_sample  title
S40     PL12  Control    GSM5102551  Control-PL12-1
S52     IL22  Control    GSM5102554  Control-IL22-1
S58     PL12  Salt       GSM5102557  Salt-PL12-1
S70     IL22  Salt       GSM5102560  Salt-IL22-1
...     (12 rows)
```

### Loading them

**Option 1: from the installed package** (no paths needed; use `"ovarian_brequinar"` for the second data set)

```python
import orthoresponse as orr
ex = orr.load_example("rice_salt")              # ex.counts, ex.samples, ex.params
res = orr.run_orthoresponse(ex.counts, ex.samples, **ex.params)
print(res)                                      # summary + asymmetry test
```

```r
library(orthoresponse)
ex <- load_example("rice_salt")                  # ex$counts, ex$samples, ex$params
res <- do.call(orthoresponse, c(list(ex$counts, ex$samples), ex$params))
res
```

**Option 2: from the files of a cloned repository** (paths relative to the repository root)

```python
import pandas as pd
import orthoresponse as orr

counts  = pd.read_csv("python/src/orthoresponse/data/rice_GSE167342_counts.tsv.gz", sep="\t", index_col=0)
samples = pd.read_csv("python/src/orthoresponse/data/rice_GSE167342_samples.tsv", sep="\t", index_col=0)

res = orr.run_orthoresponse(counts, samples, group_col="line", condition_col="treatment",
                            group_a="PL12", group_b="IL22", control="Control", stress="Salt")
```

```r
library(orthoresponse)
counts  <- as.matrix(read.delim("r-package/inst/extdata/rice_GSE167342_counts.tsv.gz",
                                row.names = 1, check.names = FALSE))
samples <- read.delim("r-package/inst/extdata/rice_GSE167342_samples.tsv", row.names = 1)

res <- orthoresponse(counts, samples, group_col = "line", condition_col = "treatment",
                     group_a = "PL12", group_b = "IL22", control = "Control", stress = "Salt")
```

For the ovarian data set, read `ovarian_GSE283381_counts.tsv.gz` and
`ovarian_GSE283381_samples.tsv` the same way and use `group_a="CisS"`,
`group_b="CisR"`, `control="Vehicle"`, `stress="Brequinar"`.

**Option 3: copy them to your own folder** and use them as a template for your data:
`orthoresponse example --out my_folder` (command line; add `--name ovarian_brequinar`
for the second set), `orr.export_example("my_folder")` (Python) or
`export_example("my_folder")` (R). This writes `my_folder/counts.tsv.gz` and
`my_folder/samples.tsv`.

The values in `line` and `treatment` are exactly what you pass as
`group_a`/`group_b` and `control`/`stress` (case-sensitive), and the sample
table must be tab-separated.

## Quick start

```python
import pandas as pd, orthoresponse as orr

counts  = pd.read_csv("counts.tsv",  sep="\t", index_col=0)   # raw counts, genes x samples
samples = pd.read_csv("samples.tsv", sep="\t", index_col=0)   # columns: line, treatment

res = orr.run_orthoresponse(counts, samples, group_col="line", condition_col="treatment",
                            group_a="sensitive", group_b="tolerant",
                            control="control", stress="drought")
print(res)           # summary table + asymmetry test (same output as R)
res.save("results")  # tables + figures (model_patterns, model_summary; PNG and PDF)
```

```r
library(orthoresponse)
res <- orthoresponse(counts, samples, group_col = "line", condition_col = "treatment",
                     group_a = "sensitive", group_b = "tolerant",
                     control = "control", stress = "drought")
res
save_orthoresponse(res, "results")   # tables + figures
```

```bash
orthoresponse run --counts counts.tsv --samples samples.tsv --group-col line \
  --condition-col treatment --group-a sensitive --group-b tolerant \
  --control control --stress drought --out results
```

Already have differential-expression results? Pass the five contrasts to
`classify_responses()`. Working with two species? Use
`reciprocal_best_hits()` and `merge_ortholog_counts()` first.

## Documentation

* **[User manual (PDF)](docs/orthoresponse_manual.pdf)**: illustrated manual with diagrams,
  figures and worked examples (LaTeX source in `docs/latex/`, build with
  `latexmk -xelatex orthoresponse_manual.tex`).
* **[User manual (Markdown)](docs/manual.md)**: concepts, input formats, decision tree,
  options, outputs, cross-species workflow, worked examples, caveats.
* R vignette: `vignette("orthoresponse")`.
* Worked examples with public, already-counted data (no alignment needed):
  * [`examples/rice_salt_GSE167342`](examples/rice_salt_GSE167342): rice,
    salt-sensitive PL12 vs salt-tolerant *Saltol* introgression line IL22.
  * [`examples/ovarian_brequinar_GSE283381`](examples/ovarian_brequinar_GSE283381):
    human ovarian cancer cells, cisplatin-resistant vs sensitive, vehicle vs brequinar.

## Repository layout

```text
python/        Python package (pip-installable, src layout, pytest suite)
r-package/     R package (roxygen2 docs, testthat suite, vignette)
docs/          user manual and figures
examples/      runnable examples in both languages
tools/         maintenance scripts (shared test fixture)
```

## Citation

If you use orthoresponse, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff) (GitHub's *Cite this repository* button).

## License

MIT, see [LICENSE](LICENSE).
