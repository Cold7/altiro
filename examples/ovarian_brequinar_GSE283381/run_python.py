"""Human ovarian cancer: OVCAR5 cisplatin-sensitive cells (group A) vs their cisplatin-resistant derivative
(group B), treated with vehicle or the DHODH inhibitor brequinar (blocks de novo pyrimidine synthesis).

Data: GEO GSE283381 (raw HTSeq counts deposited by the authors, 3 replicates), bundled with the package.
Requires: pip install "orthoresponse[all]"
Run from this folder:  python run_python.py
"""

from pathlib import Path

import orthoresponse as orr
from orthoresponse.plotting import plot_model_patterns

HERE = Path(__file__).parent


def main():
    ex = orr.load_example("ovarian_brequinar")   # counts: 59,453 genes x 12 samples
    print(ex.samples[["line", "treatment"]])

    # ex.params = group_col="line", condition_col="treatment", group_a="CisS", group_b="CisR",
    #             control="Vehicle", stress="Brequinar"
    res = orr.run_orthoresponse(ex.counts, ex.samples, **ex.params)
    out = res.save(HERE / "results_python")

    print(res.summary.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print("\nStress-ready asymmetry among single-group responders:")
    print(res.asymmetry.to_string(index=False))

    strict = orr.classify_responses(res.contrasts, preset="strict", group_names=("CisS", "CisR"))
    print("\nStress-ready genes, paper vs strict preset:",
          int((res.classification.model == "stress_ready").sum()), "->",
          int((strict.model == "stress_ready").sum()))

    # Genes the resistant cells express constitutively at the brequinar-induced level of the sensitive cells
    sr3 = strict.query("label == 'stress_ready_3'")
    top = sr3.sort_values("a_response_log2FoldChange", ascending=False).head(20)
    print("\nTop stress_ready_3 genes (strict preset):")
    print(top[["a_response_log2FoldChange", "basal_log2FoldChange"]].round(2).to_string())
    top.to_csv(out / "top_stress_ready_3_strict.tsv", sep="\t")

    fig = plot_model_patterns(res.classification, condition_labels=("Veh", "BQ"))
    fig.savefig(out / "model_patterns.png", dpi=150)
    print(f"\nResults in {out}/")


if __name__ == "__main__":
    main()
