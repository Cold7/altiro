# Examples

Each folder contains a short README, the sample table and one script per
language, writing results to `results_python/` or `results_r/`. The count
matrices are bundled with the packages (see the main README, *Example data*),
so nothing needs to be downloaded.

| folder | organism | groups (A vs B) | conditions | replicates | source |
|---|---|---|---|---|---|
| `rice_salt_GSE167342` | *Oryza sativa* | PL12 (salt-sensitive) vs IL22 (*Saltol* introgression line, tolerant) | control vs 24 h 80 mM NaCl | 3 | GEO GSE167342, `GSE167342_Raw_counts.xlsx`; `load_example("rice_salt")` |
| `ovarian_brequinar_GSE283381` | *Homo sapiens* | OVCAR5 cisplatin-sensitive vs cisplatin-resistant | vehicle vs brequinar (DHODH inhibitor) | 3 | GEO GSE283381, `GSE283381_RAW.tar` (raw HTSeq counts); `load_example("ovarian_brequinar")` |

```bash
pip install "orthoresponse[all] @ git+https://github.com/Cold7/orthoresponse.git#subdirectory=python"
cd rice_salt_GSE167342 && python run_python.py
```

```r
setwd("rice_salt_GSE167342"); source("run_r.R")
```

See section 10 of the [manual](../docs/manual.md) for the results and their
interpretation.
