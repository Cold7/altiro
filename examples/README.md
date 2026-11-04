# Examples

Each folder contains a sample sheet and one script per language. The scripts
download the public count matrix on first run (into `data/`) and write results
to `results_python/` or `results_r/`.

| folder | organism | groups (A vs B) | conditions | replicates | source |
|---|---|---|---|---|---|
| `rice_salt_GSE167342` | *Oryza sativa* | PL12 (salt-sensitive) vs IL22 (*Saltol* introgression line, tolerant) | control vs 24 h 80 mM NaCl | 3 | GEO GSE167342, `GSE167342_Raw_counts.xlsx` |
| `human_osimertinib_GSE232890` | *Homo sapiens* | PC-9 empty vector (sensitive) vs PC-9 TP53 R248Q (resistant) | DMSO vs 24 h osimertinib | 2 | GEO GSE232890, NCBI-generated raw counts |

```bash
pip install "orthoresponse[all] @ git+https://github.com/Cold7/orthoresponse.git#subdirectory=python"
cd rice_salt_GSE167342 && python run_python.py
```

```r
# install.packages("readxl")   # needed by the rice example
setwd("rice_salt_GSE167342"); source("run_r.R")
```

Download links used by the scripts:

* https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167342/suppl/GSE167342_Raw_counts.xlsx
* https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE232890&format=file&file=GSE232890_raw_counts_GRCh38.p13_NCBI.tsv.gz
* https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz

If NCBI answers with a web page instead of the file (a browser check), open the
link in a browser and save the file into the example's `data/` folder with the
name shown in the error message.

See section 10 of the [manual](../docs/manual.md) for the results and their
interpretation.
