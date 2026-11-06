# Rice salt stress (GEO GSE167342)

Leaves of the salt-sensitive *japonica* cultivar **PL12** (group A) and the
salt-tolerant introgression line **IL22** (group B; FL478 x PL12, carries the
*Saltol* QTL), control vs 24 h 80 mM NaCl, 3 replicates per group.

The data files are not in this folder: they are bundled with the packages, so
the scripts work after a plain installation.

| file | location in the repository |
|---|---|
| raw counts (36,844 genes x 12 samples) | `python/src/altiro/data/rice_GSE167342_counts.tsv.gz` (also `r-package/inst/extdata/`) |
| sample table | `python/src/altiro/data/rice_GSE167342_samples.tsv` (also `r-package/inst/extdata/`); `samples.tsv` in this folder is the same table |

Load them with `alt.load_example("rice_salt")` (Python) or
`load_example("rice_salt")` (R), or read the files directly; see the main
README, section *Example data*.

```bash
python run_python.py     # writes results_python/
Rscript run_r.R          # writes results_r/
```

Source: `GSE167342_Raw_counts.xlsx` in the GEO supplementary files
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167342). Columns S40-S72
are the authors' library names; each GEO sample (GSM) lists its library name
in the *Description* field.
