# orthoresponse user manual

Version 0.1.0 · Python and R

`orthoresponse` classifies genes (within one species) or ortholog pairs (between
two species) into **expression response models** in a 2 x 2 design: two
**groups** (species, genotypes, cultivars, cell lines...) measured under a
**control** and a **stress** condition. The default rules reproduce exactly the
decision tree of

> Contreras-Riquelme JS, Contreras M, Moyano TC, Sjoberg R, Jimenez-Gomez J,
> Alvarez JM (2025). Desert-adapted tomato *Solanum pennellii* exhibit unique
> regulatory elements and stress-ready transcriptome patterns to drought.
> *PLoS One* 20(5): e0324724. https://doi.org/10.1371/journal.pone.0324724

The Python and R implementations share the same function names, options and
outputs, and give identical classifications for identical inputs.

---

## Contents

1. [Concepts](#1-concepts)
2. [Installation](#2-installation)
3. [Quick start](#3-quick-start)
4. [Input data](#4-input-data)
5. [What the workflow does](#5-what-the-workflow-does)
6. [The decision tree](#6-the-decision-tree)
7. [Presets and options](#7-presets-and-options)
8. [Outputs](#8-outputs)
9. [Cross-species analyses (orthologs)](#9-cross-species-analyses-orthologs)
10. [Worked examples](#10-worked-examples)
11. [Reproducing the original publication](#11-reproducing-the-original-publication)
12. [Python and R: what is identical and what is not](#12-python-and-r-what-is-identical-and-what-is-not)
13. [Interpretation and statistical caveats](#13-interpretation-and-statistical-caveats)
14. [Function reference](#14-function-reference)
15. [Troubleshooting](#15-troubleshooting)
16. [Citation and license](#16-citation-and-license)

---

## 1. Concepts

### Groups and conditions

| Symbol | Meaning | Example (original publication) |
|---|---|---|
| **A** | reference group, usually the sensitive one | *Solanum lycopersicum* |
| **B** | comparison group, usually the tolerant / adapted one | *Solanum pennellii* |
| control | baseline condition | well watered |
| stress | treatment | 11 days of drought |

Which group you call A and which B only changes the case numbers, not the
biology: stress-ready cases 1-2 always describe **A primed**, cases 3-4
**B primed**.

### The five models

![Expression models](figures/model_patterns.png)

| label | description | primed group |
|---|---|---|
| `stress_ready_1` | B is induced under stress up to the control level of A | A |
| `stress_ready_2` | B is repressed under stress down to the control level of A | A |
| `stress_ready_3` | A is induced under stress up to the control level of B | B |
| `stress_ready_4` | A is repressed under stress down to the control level of B | B |
| `shared_1` / `shared_2` | induced / repressed in both groups | - |
| `unique_1` / `unique_2` | induced / repressed only in A | - |
| `unique_3` / `unique_4` | induced / repressed only in B | - |
| `opposite_1` | induced in A, repressed in B | - |
| `opposite_2` | repressed in A, induced in B | - |
| `no_response_1` / `_2` | no response; B basal level higher / lower than A | - |
| `no_response_3` | no response; no basal difference detected | - |

Two extra labels exist only with non-default options:
`unresolved_<case>` (equivalence test inconclusive, see 7.3) and
`unclassified` (missing statistics, see 7.1).

A **stress-ready** (or *primed*, *pre-adapted*) gene is one whose constitutive
expression in one group already matches the stress-induced (or stress-repressed)
expression of the other group.

---

## 2. Installation

Replace `Cold7` with the GitHub account that hosts the repository.

### Python (>= 3.9)

```bash
# everything: PyDESeq2, plotting and Excel support
pip install "orthoresponse[all] @ git+https://github.com/Cold7/orthoresponse.git#subdirectory=python"

# classifier only (numpy, pandas, scipy) - enough if you already have contrast tables
pip install "git+https://github.com/Cold7/orthoresponse.git#subdirectory=python"
```

Optional extras: `deseq` (PyDESeq2), `plot` (matplotlib), `excel` (openpyxl),
`all`, `test` (pytest). Check the installation with `orthoresponse --version`.

From a local clone (development mode):

```bash
git clone https://github.com/Cold7/orthoresponse.git
pip install -e "orthoresponse/python[all,test]"
pytest orthoresponse/python/tests
```

### R (>= 4.1)

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("DESeq2")
remotes::install_github("Cold7/orthoresponse", subdir = "r-package",
                        build_vignettes = TRUE)
library(orthoresponse)
vignette("orthoresponse")
```

From a local clone: `R CMD INSTALL orthoresponse/r-package`.

---

## 3. Quick start

### Python

```python
import pandas as pd
import orthoresponse as orr

counts = pd.read_csv("counts.tsv", sep="\t", index_col=0)    # raw counts, genes x samples
samples = pd.read_csv("samples.tsv", sep="\t", index_col=0)  # one row per sample

res = orr.run_orthoresponse(
    counts, samples,
    group_col="line", condition_col="treatment",
    group_a="PL12", group_b="IL22",          # A = sensitive, B = tolerant
    control="Control", stress="Salt",
)
print(res.summary)          # counts per model and case
print(res.asymmetry)        # is one group more often primed?
res.save("results")         # TSV files + contrast tables + normalized counts
```

### Command line (installed with the Python package)

```bash
orthoresponse run --counts counts.tsv --samples samples.tsv \
  --group-col line --condition-col treatment \
  --group-a PL12 --group-b IL22 --control Control --stress Salt \
  --out results

orthoresponse models --group-names PL12 IL22     # print the model definitions
orthoresponse classify --help                     # classify precomputed contrasts
orthoresponse rbh --help                          # reciprocal best BLAST hits
orthoresponse merge --help                        # merge two species' counts
```

### R

```r
library(orthoresponse)
counts  <- read.delim("counts.tsv", row.names = 1, check.names = FALSE)
samples <- read.delim("samples.tsv", row.names = 1)

res <- orthoresponse(counts, samples,
                     group_col = "line", condition_col = "treatment",
                     group_a = "PL12", group_b = "IL22",
                     control = "Control", stress = "Salt")
res                               # printed summary
save_orthoresponse(res, "results")
plot_model_patterns(res$classification)
```

### Figure labels

`plot_model_patterns()` labels the two conditions "C" and "S" by default.
Choose your own with `condition_labels` (control first, stress second), rename
any of the row labels with `model_labels` (keys `stress_ready`, `shared`,
`unique`, `opposite`, `no_response`; the others keep their defaults), and
optionally replace the title. Shared and opposite responses are drawn in the
first and fourth columns.

```python
from orthoresponse.plotting import plot_model_patterns
fig = plot_model_patterns(res.classification, condition_labels=("Ctrl", "NaCl"))
fig = plot_model_patterns(res.classification, condition_labels=("C", "D"),   # D = drought
                          group_names=("S. lycopersicum", "S. pennellii"), title="Tomato drought",
                          model_labels={"stress_ready": "Primed", "no_response": "Unchanged"})
fig.savefig("model_patterns.png", dpi=150)
```
```r
plot_model_patterns(res$classification, condition_labels = c("Ctrl", "NaCl"))
plot_model_patterns(res$classification, condition_labels = c("C", "D"),
                    group_names = c("S. lycopersicum", "S. pennellii"), title = "Tomato drought",
                    model_labels = c(stress_ready = "Primed", no_response = "Unchanged"))
```

---

## 4. Input data

### 4.1 Count matrix

* **Raw integer counts** (not TPM, FPKM or normalized values), features in rows,
  samples in columns, feature ids in the first column / row names.
* Estimated counts (e.g. Salmon/Kallisto `est_counts`) are accepted with
  `round_counts=True` / `round_counts = TRUE`.
* For two species, first build an ortholog-level matrix (section 9).

### 4.2 Sample table

One row per sample (row names = count-matrix column names) and at least two
columns: the group and the condition. Samples with other levels are dropped
with a warning, so one table can describe a larger experiment.

```text
sample  line  treatment
S40     PL12  Control
S52     IL22  Control
S58     PL12  Salt
S70     IL22  Salt
...
```

Every one of the four cells needs at least one sample; at least 3 replicates
per cell are strongly recommended (section 13).

### 4.3 Precomputed contrast tables (alternative entry point)

If you already ran a differential-expression analysis, call
`classify_responses()` directly with the five contrasts:

| key | numerator vs denominator |
|---|---|
| `a_response` | A stress vs A control |
| `b_response` | B stress vs B control |
| `basal` | B control vs A control |
| `b_stress_vs_a_control` | B stress vs A control |
| `a_stress_vs_b_control` | A stress vs B control |

Each table needs `log2FoldChange` and `padj` (adjusted p-value), plus `lfcSE`
if you use the equivalence test. Other column names can be mapped with
`read_contrasts(..., lfc_col=, padj_col=, se_col=)`. The sign convention
matters: a positive log2 fold change means *higher in the numerator*.

You can also pass a single wide table whose columns are
`<key>_log2FoldChange`, `<key>_padj` (and `<key>_lfcSE`).

---

## 5. What the workflow does

`run_orthoresponse()` (Python) / `orthoresponse()` (R):

1. **Design.** Samples are mapped to four cells, `A_control`, `A_stress`,
   `B_control`, `B_stress`.
2. **Filtering.** Features with a total count `<= min_total_count` (default
   10, as in the published script) are removed.
3. **Model.** One negative-binomial GLM is fitted to *all* samples with the
   design `~ cell` (DESeq2 in R, PyDESeq2 in Python). Size factors and
   dispersions are therefore estimated across both groups, as in the
   publication.
4. **Contrasts.** The five contrasts of section 4.3 are extracted with
   Benjamini-Hochberg adjustment, Cook's outlier filtering and independent
   filtering (both can be switched off). `filter_alpha` sets the level used by
   the independent filtering; by default it equals `alpha`. The published
   script used 0.01 for filtering and 0.05 for classification
   (`alpha=0.05, filter_alpha=0.01`).
5. **Classification** with the decision tree (section 6).
6. **Summaries**: counts per model with explicit denominators, and the
   primed-asymmetry test.

---

## 6. The decision tree

A contrast is **significant** when `padj < alpha` (and
`|log2FC| >= min_abs_lfc`, default 0). Its **direction** is the sign of the
log2 fold change rounded to `sign_digits` (default 2) decimals.

```text
                    A responds? (a_response significant)
               no  /                                 \ yes
   B responds? (b_response)                    B responds? (b_response)
      no /          \ yes                        no /            \ yes
 NO RESPONSE      B up or down?             A up or down?     SHARED / OPPOSITE
 basal signif.?   B stress vs A control     A stress vs B control   by the two signs:
 +  -> NR 1       not signif.  -> SR 1 (up) not signif. -> SR 3 (up)   + +  -> shared 1
 -  -> NR 2                       SR 2 (down)              SR 4 (down) - -  -> shared 2
 ns -> NR 3       significant  -> UR 3 (up) significant -> UR 1 (up)   + -  -> opposite 1
                                  UR 4 (down)              UR 2 (down) other -> opposite 2
```

(SR = stress-ready, UR = unique response, NR = no response.)

### Flags

The `flags` column records why a call may deserve attention:

* `na_treated_as_not_de`: at least one contrast consulted for this feature had
  a missing adjusted p-value and, under `na_action="not_de"`, it was read as
  "not significant". Missing values come from DESeq2's independent filtering
  (low counts) and Cook's outlier detection.
* `basal_difference_not_significant`: the feature would be stress-ready under
  the paper rules but was demoted to unique because no significant basal
  difference supports it (only with `require_basal_difference`).

---

## 7. Presets and options

| option | `paper` (default) | `strict` | meaning |
|---|---|---|---|
| `alpha` | 0.05 | 0.05 | adjusted p-value threshold |
| `min_abs_lfc` | 0 | 0 | minimum absolute log2FC for "significant" |
| `na_action` | `not_de` | `unclassified` | how missing `padj` is handled |
| `require_basal_difference` | False | True | stress-ready needs a basal difference |
| `equivalence_margin` | none | none | TOST equivalence for "same level" |
| `sign_digits` | 2 | 2 | rounding before taking the sign |

Any option can be combined with any preset:

```python
orr.classify_responses(tables, preset="strict", equivalence_margin=0.5)
```
```r
classify_responses(tables, preset = "strict", equivalence_margin = 0.5)
```

### 7.1 Missing adjusted p-values (`na_action`)

Under the paper rules a missing value counts as "not different". In the
contrasts that test whether two levels are *equal*, this pushes features
towards stress-ready. With `na_action="unclassified"` such features get the
label `unclassified` instead. Most of them are low-count genes removed by the
independent filtering of one contrast.

### 7.2 Basal difference (`require_basal_difference`)

If A is induced up to B's control level (stress-ready 3), B's control level must
be higher than A's control level. The paper rules never check this: they only
check that A-stress is *not significantly different* from B-control. When the
two groups have the same basal level, a gene that responds in A alone (a unique
response) can pass that check simply because the induction is modest. With
this option, a stress-ready call also requires a significant `basal` contrast in
the expected direction (negative for cases 1 and 4, positive for cases 2 and 3);
otherwise the feature becomes the corresponding unique response.

### 7.3 Equivalence (`equivalence_margin`)

"Not significantly different" is not evidence of equality. With
`equivalence_margin = m` (log2 units), the matching contrast must pass a two
one-sided test (TOST) for `|log2FC| < m`:

`p = max( P(Z > (lfc + m)/se), P(Z < (lfc - m)/se) )`, BH-adjusted across
features, which is the construction DESeq2 uses with
`altHypothesis = "lessAbs"`. Features that are neither different nor
equivalent are labelled `unresolved_<case>`. Choose `m` before looking at the
results (0.5 to 1 is common); with 2-3 replicates few genes can be shown to be
equivalent within 0.5 log2 units, which is itself informative.

### 7.4 Which preset should I use?

* To reproduce or compare with the original publication: `paper`.
* To claim that a group is pre-adapted, or when the two groups are
  near-isogenic (e.g. an introgression line and its recurrent parent): `strict`,
  ideally with an equivalence margin, and report both.

---

## 8. Outputs

### 8.1 `classification`

One row per feature (index / row names = feature id):

| column | content |
|---|---|
| `label` | e.g. `stress_ready_3` |
| `model`, `case` | the two parts of the label |
| `primed_group` | group whose control level equals the other's stressed level (stress-ready only) |
| `responding_group` | A, B, `both` or `none` |
| `description` | human-readable sentence with your group names |
| `flags` | see section 6 |
| `<key>_log2FoldChange`, `<key>_padj` | statistics of the five contrasts |

### 8.2 `summary`

`label, model, case, n, pct_of_responsive, pct_of_all`. Two denominators are
reported on purpose:

* `pct_of_responsive`: features that respond in at least one group
  (stress-ready, shared, unique, opposite, unresolved). The percentages in
  Fig. 3 of the original publication use this denominator.
* `pct_of_all`: every classified feature.

### 8.3 `asymmetry` (`primed_asymmetry()`)

Among features that respond in only one group, the fraction that ends in a
stress-ready state, for A-only and B-only responders, compared with Fisher's
exact test. A larger raw number of stress-ready genes in one direction can
simply reflect that one group has more responsive genes; this test removes that
effect. The odds ratio reported is the sample odds ratio
`(SR_A * U_B) / (U_A * SR_B)`.

### 8.4 Files written by `save()` / `save_orthoresponse()`

```text
results/
  classification.tsv
  summary.tsv
  primed_asymmetry.tsv
  normalized_counts.tsv
  contrasts/contrast_<key>.tsv     (full DESeq2 / PyDESeq2 tables)
```

---

## 9. Cross-species analyses (orthologs)

1. **Quantify each species against its own reference** (any aligner or
   pseudo-aligner) to obtain raw gene-level counts.
2. **Define one-to-one orthologs.** The package implements reciprocal best
   BLASTp hits, the method of the original publication:

   ```bash
   makeblastdb -in B.faa -dbtype prot -out B
   makeblastdb -in A.faa -dbtype prot -out A
   blastp -query A.faa -db B -outfmt 6 -evalue 1e-5 -max_target_seqs 5 -num_threads 8 > A_vs_B.tsv
   blastp -query B.faa -db A -outfmt 6 -evalue 1e-5 -max_target_seqs 5 -num_threads 8 > B_vs_A.tsv
   ```

   ```python
   orth = orr.reciprocal_best_hits("A_vs_B.tsv", "B_vs_A.tsv", strip_pattern=r"\.\d+$")
   ```
   ```r
   orth <- reciprocal_best_hits("A_vs_B.tsv", "B_vs_A.tsv", strip_pattern = "\\.\\d+$")
   ```

   `strip_pattern` removes isoform / version suffixes so that protein ids match
   the gene ids of the count matrices. The best hit is chosen by bitscore, then
   e-value, then identifier. Any other one-to-one ortholog table (OrthoFinder
   single-copy orthogroups, Ensembl Compara, PLAZA...) works as long as it has
   one column per species.
3. **Merge the two count matrices:**

   ```python
   counts = orr.merge_ortholog_counts(counts_a, counts_b, orth)
   ```

   Rows become ortholog ids; columns are the samples of A followed by B
   (sample names must be unique). `length_normalize=True` applies the
   correction of the original publication, `Q = round(E / L_kb)`, using
   per-gene lengths. It puts orthologs of different lengths on one scale but
   changes the count distribution DESeq2 assumes, so use it to replicate that
   analysis or as a sensitivity check rather than by default.
4. **Run the workflow** with the species as groups.

Basal differences between species mix biology with technical effects
(reference quality, transcript length, mappability, annotation of UTRs). Treat
`no_response_1/2` and the basal part of stress-ready calls with that in mind.

---

## 10. Worked examples

Both examples start from public, already counted matrices; no read alignment is
needed. Scripts are in `examples/` (one Python and one R script each); they
download the data on first use. Numbers below come from the R/DESeq2 run.

### 10.1 Plant: rice salt stress (GEO GSE167342)

* **Design.** Leaves of the salt-sensitive *japonica* cultivar **PL12** (A) and
  the salt-tolerant introgression line **IL22** (B; derived from FL478 x PL12,
  carries the *Saltol* QTL), control vs 24 h 80 mM NaCl, 3 replicates.
* **Data.** `GSE167342_Raw_counts.xlsx` from the GEO supplementary files
  (gene-level raw counts, 36,844 genes x 12 samples; columns S40-S72 are mapped
  to samples in `examples/rice_salt_GSE167342/samples.tsv`).
  Publication: PMID 35126422.

```bash
cd examples/rice_salt_GSE167342
python run_python.py     # or: Rscript run_r.R
```

Results (25,045 genes after filtering):

| model | paper | strict |
|---|---|---|
| stress-ready 1 / 2 / 3 / 4 | 281 / 407 / 1,744 / 1,793 | 35 / 72 / 297 / 189 |
| shared 1 / 2 | 563 / 335 | 563 / 335 |
| unique 1 / 2 / 3 / 4 | 617 / 628 / 152 / 185 | 2,024 / 2,180 / 389 / 514 |
| opposite 1 / 2 | 88 / 41 | 88 / 41 |
| no response 1 / 2 / 3 | 174 / 196 / 17,841 | 174 / 196 / 13,587 |
| unclassified | 0 | 4,361 |

![Rice model patterns](figures/rice_GSE167342_model_patterns.png)

Lessons from this data set:

* PL12 (sensitive) changes many more genes than IL22 (4,782 vs 1,025
  single-line responses), consistent with IL22 experiencing less stress.
* IL22 and PL12 are nearly isogenic, so basal differences are rare (only 370
  non-responsive genes differ between the lines).
  Under the paper rules, 62 % of responsive genes are called stress-ready,
  mostly because a gene that responds in one line is "not significantly
  different" from the unchanged basal level of the other. Requiring a basal
  difference (`strict`) keeps 593 of 4,225.
* The paper-rule asymmetry (74.0 % of PL12-only responses end at IL22's level vs
  67.1 % in the other direction, p = 1.2e-5) disappears under `strict`
  (10.4 % vs 10.6 %, p = 0.82).
* With an equivalence margin of 0.5 log2 units and 3 replicates, only one gene
  can be shown to be stress-ready; the rest become `unresolved`.

### 10.2 Human: drug tolerance in lung cancer cells (GEO GSE232890)

* **Design.** PC-9 non-small-cell lung cancer cells (EGFR mutant) expressing an
  empty vector (**EV**, A, osimertinib-sensitive) or the TP53 gain-of-function
  mutant **R248Q** (B, early osimertinib resistance), treated 24 h with DMSO
  (control) or osimertinib (stress), 2 replicates.
* **Data.** Raw counts generated by NCBI from the submitted reads
  (`GSE232890_raw_counts_GRCh38.p13_NCBI.tsv.gz`, NCBI Gene ids). Gene symbols
  come from NCBI's `Homo_sapiens.gene_info.gz`. NCBI produces such matrices for
  most human (and mouse) RNA-seq series in GEO: on any series page, open
  "Download RNA-seq counts", or search GEO DataSets with the filter
  `"rnaseq counts"[Filter]`.

```bash
cd examples/human_osimertinib_GSE232890
python run_python.py     # or: Rscript run_r.R
```

Results (21,583 genes): most responsive genes are shared (6,156 of 8,794),
1,049 are stress-ready under the paper rules and 226 under `strict`. The
asymmetry test does not support a globally primed resistant line
(paper: 45.0 % of R248Q-only vs 37.2 % of EV-only responses end at the other
line's level; strict: 9.1 % vs 8.5 %, p = 0.61). Stress-ready 3 genes (induced
by osimertinib in EV and constitutively high in R248Q) are listed in
`results_*/top_stress_ready_3.tsv`. With only 2 replicates per group, treat
gene-level calls as exploratory.

### 10.3 Finding other data sets

The method needs a full 2 x 2 design: both groups under both conditions. Good
sources of public, already-counted matrices:

* **Human / mouse**: NCBI-generated counts in GEO (see above). Typical designs:
  resistant vs parental cell line x vehicle vs drug; knockout vs wild type x
  treatment; two tissues or donors x stimulus.
* **Plants and other organisms**: GEO supplementary files named
  `*count*`; search GEO DataSets for "tolerant sensitive drought", "salt
  tolerant genotype", etc., and check that raw counts (not FPKM/TPM) are given.
  Other examples found while preparing this manual: sorghum drought,
  two contrasting cultivars (GSE276121, raw counts with gene lengths).

---

## 11. Reproducing the original publication

The published script wrote five CSV files with DESeq2 results. Classify them
directly:

```python
paths = {
    "a_response":            "Sly.drought_Sly.control.csv",
    "b_response":            "Spenn.drought_Spenn.control.csv",
    "basal":                 "Spenn.control_Sly.control.csv",
    "b_stress_vs_a_control": "Spenn.drought_Sly.control.csv",
    "a_stress_vs_b_control": "Sly.drought_Spenn.control.csv",
}
tabs = orr.read_contrasts(paths, id_col="Row.names")
cls = orr.classify_responses(tabs, group_names=("S. lycopersicum", "S. pennellii"))
orr.summarize_classification(cls)
```

The Python test suite contains a verbatim copy of the published notebook logic
(`python/tests/reference_notebook.py`) and checks that `preset="paper"` gives
the same label for every feature on randomly generated tables that cover every
branch, including missing values and zero fold changes.

To refit from counts with the published settings use
`min_total_count=10, alpha=0.05, filter_alpha=0.01`, and
`merge_ortholog_counts(..., length_normalize=True)` for the length correction.

---

## 12. Python and R: what is identical and what is not

* **Classifier**: identical. Given the same contrast tables, both languages
  return the same label for every feature (tested on 240,000 random
  classifications with four option sets, on a shared fixture in both test
  suites, and on the two worked examples).
* **Differential expression engine**: DESeq2 (R) and PyDESeq2 (Python) are
  separate implementations. Fold changes agree almost perfectly
  (correlation > 0.9999 on the examples), but adjusted p-values near the
  threshold can differ, so a few borderline genes change label: 99.5 % of
  labels agreed on the rice example (3 replicates) and 96.1 % on the human one
  (2 replicates). For exact cross-language agreement, compute the contrasts
  once (e.g. in R), write them with `write_contrasts()`, and classify them in
  either language.

---

## 13. Interpretation and statistical caveats

* **Absence of evidence.** Under the paper rules, "same level" means "not
  significantly different". Low power (few replicates, low counts) therefore
  inflates stress-ready calls. Report the `strict` preset alongside, or use an
  equivalence margin.
* **False positives flow into stress-ready.** A false-positive response in one
  group, for a gene whose basal level does not differ between groups, is
  classified as stress-ready under the paper rules. `require_basal_difference`
  blocks this route.
* **Unique vs shared.** "Significant in A and not in B" does not prove that the
  responses differ. If the difference of responses itself is your question,
  test the interaction (`~ group + group:condition`) separately.
* **Directional claims.** Use `primed_asymmetry()` before concluding that one
  group is "more pre-adapted": raw counts of stress-ready cases depend on how
  many genes respond in each group.
* **Equal stress?** The same treatment may be perceived differently by the two
  groups (e.g. a tolerant genotype uses less water). Fewer responses in one
  group can mean less perceived stress rather than pre-adaptation.
* **Replicates.** Use at least 3 per cell; with 2, treat results as
  exploratory.

---

## 14. Function reference

| Python | R | purpose |
|---|---|---|
| `run_orthoresponse()` | `orthoresponse()` | counts -> contrasts -> classification |
| `run_contrasts()` | `run_contrasts()` | fit (Py)DESeq2 and extract the five contrasts |
| `build_design()` | `build_design()` | map samples to the four cells |
| `classify_responses()` | `classify_responses()` | decision tree |
| `classification_options()` | `classification_options()` | presets and options |
| `model_definitions()` | `model_definitions()` | table of models and cases |
| `summarize_classification()` | `summarize_classification()` | counts and percentages |
| `primed_asymmetry()` | `primed_asymmetry()` | direction test for stress-ready |
| `equivalence_test()` | `equivalence_test()` | TOST for log2 fold changes |
| `read_contrasts()` / `write_contrasts()` | same | contrast tables from / to files |
| `reciprocal_best_hits()` | `reciprocal_best_hits()` | orthologs from two BLAST files |
| `read_blast_tabular()`, `best_hits()` | same | BLAST helpers |
| `merge_ortholog_counts()` | `merge_ortholog_counts()` | ortholog-level count matrix |
| `plotting.plot_model_patterns()` | `plot_model_patterns()` | the 15 patterns with counts; custom `condition_labels`, `model_labels`, `group_names`, `title` |
| `plotting.plot_summary()` | `plot_summary()` | bar chart of models |
| `plotting.plot_feature()` | - | normalized counts of one feature |
| `OrthoResponseResult.save()` | `save_orthoresponse()` | write all outputs |
| `CONTRAST_KEYS` | `contrast_keys()` | the five contrast names |

Every function is documented in place: `help(orr.classify_responses)` in
Python, `?classify_responses` in R.

---

## 15. Troubleshooting

| message | cause and fix |
|---|---|
| `Counts are not integers` | you passed TPM/FPKM or estimated counts; use raw counts or `round_counts=True` |
| `Empty design cell(s)` | a group x condition combination has no samples; check the level names |
| `Samples missing from the count matrix` | sample table row names do not match count columns |
| `Contrast tables cover different features` | the five tables were filtered differently; only shared features are used |
| `equivalence_margin requires an 'lfcSE' column` | keep the standard-error column when exporting contrasts |
| `run_contrasts() needs PyDESeq2` | `pip install "orthoresponse[deseq]"` |
| NCBI download returns HTML | NCBI sometimes asks for a browser check; download the file in a browser and place it where the script expects it |
| `The dispersion trend curve fitting did not converge` (PyDESeq2) | informative warning; PyDESeq2 falls back to a mean dispersion trend |

---

## 16. Citation and license

If you use this package, please cite the software (see `CITATION.cff`).
DESeq2 / PyDESeq2 should be cited as well:

* Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and
  dispersion for RNA-seq data with DESeq2. *Genome Biology* 15: 550.
* Muzellec B, Telenczuk M, Cabeli V, Andreux M (2023). PyDESeq2: a python
  package for bulk RNA-seq differential expression analysis.
  *Bioinformatics* 39: btad547.

License: MIT.
