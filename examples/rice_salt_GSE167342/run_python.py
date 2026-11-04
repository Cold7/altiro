"""Rice salt stress: salt-sensitive PL12 (group A) vs salt-tolerant Saltol introgression line IL22 (group B).

Data: GEO GSE167342 (leaves, control vs 24 h 80 mM NaCl, 3 replicates), bundled with the package.
Requires: pip install "orthoresponse[all]"
Run from this folder:  python run_python.py
"""

from pathlib import Path

import orthoresponse as orr

HERE = Path(__file__).parent


def main():
    ex = orr.load_example("rice_salt")        # counts: 36,844 genes x 12 samples (S40 ... S72)
    print(ex.samples[["line", "treatment"]])

    # ex.params = group_col="line", condition_col="treatment", group_a="PL12", group_b="IL22",
    #             control="Control", stress="Salt"
    res = orr.run_orthoresponse(ex.counts, ex.samples, **ex.params)
    # tables + figures (model_patterns.png/pdf, model_summary.png/pdf)
    out = res.save(HERE / "results_python", condition_labels=("Ctrl", "NaCl"))

    print(res.summary.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print("\nStress-ready asymmetry among single-group responders:")
    print(res.asymmetry.to_string(index=False))

    # Same data, stricter rules: NA -> unclassified, basal difference required
    strict = orr.classify_responses(res.contrasts, preset="strict", group_names=("PL12", "IL22"))
    print("\nStress-ready genes, paper vs strict preset:",
          int((res.classification.model == "stress_ready").sum()), "->",
          int((strict.model == "stress_ready").sum()))

    # Genes for which the tolerant line is primed (IL22 constitutively at the salt-induced level of PL12)
    primed = res.classification.query("label == 'stress_ready_3'")
    top = primed.sort_values("a_response_log2FoldChange", ascending=False).head(15)
    top[["a_response_log2FoldChange", "a_response_padj", "basal_log2FoldChange"]].to_csv(
        out / "top_stress_ready_3.tsv", sep="\t")

    print(f"\nResults in {out}/")


if __name__ == "__main__":
    main()
