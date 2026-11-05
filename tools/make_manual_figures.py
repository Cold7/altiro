"""Regenerate the figures of the LaTeX manual (docs/latex/figures).

Run the R examples first (examples/*/run_r.R), then from the repository root:
    python tools/make_manual_figures.py
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import orthoresponse as orr
from orthoresponse.plotting import plot_model_patterns, plot_summary, COLOR_A, COLOR_B

plt.rcParams.update({"font.family": "Inter", "pdf.fonttype": 42})
OUT = "docs/latex/figures"
EX = {"rice": ("examples/rice_salt_GSE167342", ("PL12", "IL22"), ("Ctrl", "NaCl"), "rice_salt"),
      "ovarian": ("examples/ovarian_brequinar_GSE283381", ("CisS", "CisR"), ("Veh", "BQ"), "ovarian_brequinar")}

plot_model_patterns(group_names=("Group A (reference)", "Group B (comparison)")).savefig(f"{OUT}/patterns_generic.pdf", bbox_inches="tight")

cls_all, norm_all = {}, {}
for key, (d, names, cond, exname) in EX.items():
    tabs = orr.read_contrasts({k: f"{d}/results_r/contrasts/contrast_{k}.tsv" for k in orr.CONTRAST_KEYS})
    cls = orr.classify_responses(tabs, group_names=names)
    strict = orr.classify_responses(tabs, preset="strict", group_names=names)
    cls_all[key] = (cls, strict)
    plot_model_patterns(cls, condition_labels=cond).savefig(f"{OUT}/{key}_patterns.pdf", bbox_inches="tight")
    f = plot_summary(cls); f.set_size_inches(5.2, 4.6); f.tight_layout(); f.savefig(f"{OUT}/{key}_summary.pdf", bbox_inches="tight")
    norm = pd.read_csv(f"{d}/results_r/normalized_counts.tsv", sep="\t", index_col=0)
    norm.index = norm.index.astype(str); norm_all[key] = norm
    plt.close("all")

# real genes, one per model (ovarian, strict preset for primed)
def gene_panel(key, picks, fname):
    d, names, cond, exname = EX[key]
    ex = orr.load_example(exname)
    design = orr.build_design(ex.samples, **{k: ex.params[k] for k in ("group_col","condition_col","group_a","group_b","control","stress")})
    norm = norm_all[key]
    fig, axes = plt.subplots(1, len(picks), figsize=(2.55 * len(picks), 2.7), sharey=False)
    cells = ["A_control", "A_stress", "B_control", "B_stress"]
    xt = [f"{names[0]}\n{cond[0]}", f"{names[0]}\n{cond[1]}", f"{names[1]}\n{cond[0]}", f"{names[1]}\n{cond[1]}"]
    for ax, (gene, label) in zip(axes, picks):
        v = np.log2(norm.loc[gene, design.index] + 1)
        for i, c in enumerate(cells):
            vals = v[design["cell"] == c]
            col = COLOR_A if c.startswith("A") else COLOR_B
            ax.scatter(np.full(len(vals), i) + np.linspace(-0.08, 0.08, len(vals)), vals, color=col, s=22, zorder=3)
            ax.hlines(vals.mean(), i - 0.28, i + 0.28, color="k", lw=1.4)
        ax.set_xticks(range(4)); ax.set_xticklabels(xt, fontsize=6.5)
        ax.set_title(f"{gene}\n{label}", fontsize=8.5)
        ax.tick_params(axis="y", labelsize=7); ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("log2(normalized count + 1)", fontsize=7.5)
    fig.tight_layout(); fig.savefig(f"{OUT}/{fname}", bbox_inches="tight"); plt.close(fig)

cls, strict = cls_all["ovarian"]
def top(df, label, col="a_response_log2FoldChange", asc=False, coding=True):
    sub = df[df.label == label]
    if coding:
        sub = sub[~sub.index.str.match(r"^(A[CLP]\d|LINC|.*-AS\d|.*-DT|SNOR|MIR|RNU|Z\d|U\d)") & ~sub.index.str.contains(r"P\d+$")]
    sub = sub[(sub.a_response_padj < 0.001) | (sub.b_response_padj < 0.001)] if label not in ("no_response_3",) else sub
    expr = norm_all["ovarian"].mean(axis=1)
    sub = sub[expr.reindex(sub.index).fillna(0) >= 300]
    return sub.sort_values(col, ascending=asc).index[0]
picks = [(top(strict, "primed_3"), "primed_3"),
         (top(cls, "shared_1"), "shared_1"),
         (top(cls, "unique_1"), "unique_1"),
         (top(cls, "unique_3", col="b_response_log2FoldChange"), "unique_3")]
print("ovarian picks:", picks)
gene_panel("ovarian", picks, "ovarian_genes.pdf")

# paper vs strict, both data sets
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.0))
for ax, key, title in zip(axes, ["rice", "ovarian"], ["Rice (GSE167342)", "Ovarian cells (GSE283381)"]):
    c, s = cls_all[key]
    labs = [f"primed_{i}" for i in range(1, 5)]
    pv = [int((c.label == l).sum()) for l in labs]; sv = [int((s.label == l).sum()) for l in labs]
    x = np.arange(4)
    ax.bar(x - 0.2, pv, 0.4, color="#9aa7c7", label="paper preset")
    ax.bar(x + 0.2, sv, 0.4, color="#2f4b8a", label="strict preset")
    for xi, a, b in zip(x, pv, sv):
        ax.text(xi - 0.2, a, str(a), ha="center", va="bottom", fontsize=6.5)
        ax.text(xi + 0.2, b, str(b), ha="center", va="bottom", fontsize=6.5)
    ax.set_xticks(x); ax.set_xticklabels([f"primed {i}" for i in range(1, 5)], fontsize=8)
    ax.set_title(title, fontsize=9); ax.spines[["top", "right"]].set_visible(False); ax.tick_params(axis="y", labelsize=7)
axes[0].set_ylabel("Primed genes", fontsize=8); axes[1].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/paper_vs_strict.pdf", bbox_inches="tight"); plt.close(fig)

# asymmetry fractions
fig, axes = plt.subplots(1, 2, figsize=(8.6, 2.8))
for ax, key, title in zip(axes, ["rice", "ovarian"], ["Rice (GSE167342)", "Ovarian cells (GSE283381)"]):
    c, s = cls_all[key]
    a_p, a_s = orr.primed_asymmetry(c), orr.primed_asymmetry(s)
    x = np.arange(2)
    ax.bar(x - 0.2, 100 * a_p.fraction_primed, 0.4, color="#9aa7c7", label="paper")
    ax.bar(x + 0.2, 100 * a_s.fraction_primed, 0.4, color="#2f4b8a", label="strict")
    ax.set_xticks(x); ax.set_xticklabels([f"{g}-only\nresponders" for g in a_p.responding_group], fontsize=7.5)
    ax.set_ylim(0, 100); ax.set_ylabel("% ending primed", fontsize=7.5); ax.tick_params(axis="y", labelsize=7)
    ax.set_title(f"{title}\nFisher p: paper {a_p.fisher_pvalue[0]:.1g}, strict {a_s.fisher_pvalue[0]:.2f}", fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)
axes[1].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/asymmetry.pdf", bbox_inches="tight"); plt.close(fig)

for key in EX:
    c, s = cls_all[key]
    print(key, len(c), orr.summarize_classification(c).set_index("label")["n"].to_dict())
