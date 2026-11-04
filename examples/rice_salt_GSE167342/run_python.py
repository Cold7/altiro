"""Rice salt stress: salt-sensitive PL12 (group A) vs salt-tolerant Saltol introgression line IL22 (group B).

Data: GEO GSE167342 (leaves, control vs 24 h 80 mM NaCl, 3 replicates).
Requires: pip install "orthoresponse[all]"
Run from this folder:  python run_python.py
"""

import urllib.request
from pathlib import Path

import pandas as pd

import orthoresponse as orr
from orthoresponse.plotting import plot_model_patterns

HERE = Path(__file__).parent
DATA = HERE / "data"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167342/suppl/GSE167342_Raw_counts.xlsx"


def download():
    DATA.mkdir(exist_ok=True)
    dest = DATA / "GSE167342_Raw_counts.xlsx"
    if not dest.exists():
        print("Downloading", URL)
        urllib.request.urlretrieve(URL, dest)
    return dest


def main():
    counts = pd.read_excel(download(), index_col=0)          # 36,844 genes x 12 samples (S40 ... S72)
    samples = pd.read_csv(HERE / "samples.tsv", sep="\t", index_col=0)

    res = orr.run_orthoresponse(
        counts, samples,
        group_col="line", condition_col="treatment",
        group_a="PL12", group_b="IL22",          # A = sensitive, B = tolerant
        control="Control", stress="Salt",
    )
    out = res.save(HERE / "results_python")

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

    fig = plot_model_patterns(res.classification)
    fig.savefig(out / "model_patterns.png", dpi=150)
    print(f"\nResults in {out}/")


if __name__ == "__main__":
    main()
