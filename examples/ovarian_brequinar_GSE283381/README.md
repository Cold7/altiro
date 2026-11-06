# Human ovarian cancer cells, cisplatin resistance x brequinar (GEO GSE283381)

OVCAR5 high-grade serous ovarian cancer cells: the cisplatin-sensitive parental
line (**CisS**, group A) and a cisplatin-resistant derivative obtained by
repeated cisplatin treatment (**CisR**, group B), each treated with vehicle
(control) or **brequinar**, an inhibitor of DHODH and therefore of de novo
pyrimidine synthesis (stress), 3 replicates per group. The authors found
pyrimidine metabolism enriched in cisplatin-resistant cells and used this
experiment to test that dependency.

The data files are bundled with the packages (not in this folder):

| file | location in the repository |
|---|---|
| raw counts (59,453 genes x 12 samples) | `python/src/altiro/data/ovarian_GSE283381_counts.tsv.gz` (also `r-package/inst/extdata/`) |
| sample table | `python/src/altiro/data/ovarian_GSE283381_samples.tsv` (also `r-package/inst/extdata/`); `samples.tsv` in this folder has the same first three columns |

Load them with `alt.load_example("ovarian_brequinar")` (Python) or
`load_example("ovarian_brequinar")` (R), or read the files directly as shown
in the main README.

```bash
python run_python.py     # writes results_python/
Rscript run_r.R          # writes results_r/
```

Source: `GSE283381_RAW.tar` in the GEO supplementary files
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283381): one raw
HTSeq-count file per sample (reads aligned with STAR to GRCh38, Ensembl
release 91), merged into one matrix without any transformation. Column names
are the authors' library names (`OV5WTC1` = sensitive, control, replicate 1;
`OV5CRBQ3` = resistant, brequinar, replicate 3; ...). Publication: PMID 41091828.
