# altiro (Python)

ALTIRO, Analysis of Latent Transcriptional Induction and Response Outcomes.

Classify genes or ortholog pairs into primed (stress-ready), shared, unique, opposite and
no-response expression models in a 2 x 2 design (two groups x control/stress).
The default rules reproduce the decision tree of Contreras-Riquelme *et al.*
(2025) *PLoS One* 20(5): e0324724.

```bash
pip install "altiro[all] @ git+https://github.com/Cold7/altiro.git#subdirectory=python"
```

```python
import altiro as alt

res = alt.run_altiro(counts, samples, group_col="line", condition_col="treatment",
                     group_a="sensitive", group_b="tolerant",
                     control="control", stress="drought")
print(res.summary)
res.save("results")
```

Command line: `altiro --help`.

Full documentation: see `docs/manual.md` in the repository
(https://github.com/Cold7/altiro/blob/main/docs/manual.md).

Extras: `deseq` (PyDESeq2), `plot` (matplotlib), `excel` (openpyxl), `all`, `test`.
