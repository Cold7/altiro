"""Human ovarian cancer: OVCAR5 cisplatin-sensitive cells (group A) vs their cisplatin-resistant derivative
(group B), treated with vehicle or the DHODH inhibitor brequinar (blocks de novo pyrimidine synthesis).

Data: GEO GSE283381 (raw HTSeq counts deposited by the authors, 3 replicates), bundled with the package.
Requires: pip install "altiro[all]"
Run from this folder:  python run_python.py
"""

from pathlib import Path

import altiro as alt

HERE = Path(__file__).parent


def main():
    ex = alt.load_example("ovarian_brequinar")   # counts: 59,453 genes x 12 samples
    print(ex.samples[["line", "treatment"]])

    # ex.params = group_col="line", condition_col="treatment", group_a="CisS", group_b="CisR",
    #             control="Vehicle", stress="Brequinar"
    res = alt.run_altiro(ex.counts, ex.samples, **ex.params)
    # tables + figures (model_patterns.png/pdf, model_summary.png/pdf)
    out = res.save(HERE / "results_python", condition_labels=("Veh", "BQ"))

    print(res.summary.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print("\nPrimed-state asymmetry among single-group responders:")
    print(res.asymmetry.to_string(index=False))

    strict = alt.classify_responses(res.contrasts, preset="strict", group_names=("CisS", "CisR"))
    print("\nPrimed genes, paper vs strict preset:",
          int((res.classification.model == "primed").sum()), "->",
          int((strict.model == "primed").sum()))

    # Genes the resistant cells express constitutively at the brequinar-induced level of the sensitive cells
    sr3 = strict.query("label == 'primed_3'")
    top = sr3.sort_values("a_response_log2FoldChange", ascending=False).head(20)
    print("\nTop primed_3 genes (strict preset):")
    print(top[["a_response_log2FoldChange", "basal_log2FoldChange"]].round(2).to_string())
    top.to_csv(out / "top_primed_3_strict.tsv", sep="\t")

    print(f"\nResults in {out}/")


if __name__ == "__main__":
    main()
