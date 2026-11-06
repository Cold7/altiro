"""Annotated figure for the 'Customizing names and labels' section of the manual.

Run the R rice example first (examples/rice_salt_GSE167342/run_r.R), then from the repository root:
    python tools/make_customization_figure.py
"""
import warnings; warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import altiro as alt
from altiro.plotting import plot_model_patterns

plt.rcParams.update({"font.family": "Inter", "pdf.fonttype": 42})
D = "examples/rice_salt_GSE167342/results_r/contrasts"
tabs = alt.read_contrasts({k: f"{D}/contrast_{k}.tsv" for k in alt.CONTRAST_KEYS})
cls = alt.classify_responses(tabs, group_names=("PL12", "IL22"))

fig = plot_model_patterns(cls, group_names=("PL12 (sensitive)", "IL22 (tolerant)"),
                          condition_labels=("Ctrl", "NaCl"),
                          model_labels={"primed": "Pre-adapted", "no_response": "Unchanged"},
                          title="Rice under salt stress: customized labels", figsize=(11, 10))
fig.subplots_adjust(left=0.25, right=0.98, top=0.90, bottom=0.05)
fig.canvas.draw()
axes = [a for a in fig.axes]
grid = [axes[i * 4:(i + 1) * 4] for i in range(5)]
inv = fig.transFigure.inverted()

def fig_xy(ax, x, y):  # axes fraction -> figure fraction
    return inv.transform(ax.transAxes.transform((x, y)))

ORANGE, NAVY = "#C76E00", "#1F3A5F"
def callout(text, xy_text, xy_target, ha="left"):
    fig.text(*xy_text, text, fontsize=9.5, family="DejaVu Sans Mono", color=NAVY, ha=ha, va="center",
             bbox=dict(boxstyle="round,pad=0.35", fc="#FFF4E5", ec=ORANGE, lw=1.1))
    fig.add_artist(matplotlib.patches.FancyArrowPatch(xy_text, xy_target, transform=fig.transFigure,
                   arrowstyle="-|>", mutation_scale=12, color=ORANGE, lw=1.3,
                   shrinkA={"left": 2, "right": 2}.get(ha, 2) + 70 if False else 4, shrinkB=3))

# title
t = fig._suptitle
fig._suptitle.set_y(0.955); fig.canvas.draw()
bb = inv.transform(t.get_window_extent())
callout("title", (0.10, 0.955), (bb[0, 0] - 0.005, 0.955))
# row labels
for r, txt in ((0, 'model_labels\n  ["primed"]'), (4, 'model_labels\n  ["no_response"]')):
    lab = grid[r][0].yaxis.label
    lb = inv.transform(lab.get_window_extent())
    callout(txt, (0.02, (lb[0, 1] + lb[1, 1]) / 2), (lb[0, 0] - 0.003, (lb[0, 1] + lb[1, 1]) / 2))
# condition ticks
tick = grid[1][0].get_xticklabels()[0]
tb = inv.transform(tick.get_window_extent())
callout("condition_labels", (0.02, tb[1, 1] + 0.035), (tb[0, 0], (tb[0, 1] + tb[1, 1]) / 2))
# legend
leg = fig.legends[0]
lgb = inv.transform(leg.get_window_extent())
callout("group_names", (lgb[0, 0] - 0.02, lgb[1, 1] + 0.06), ((lgb[0, 0] + lgb[1, 0]) / 2, lgb[1, 1]), ha="right")
# counts
xl = grid[3][3].xaxis.label
xb = inv.transform(xl.get_window_extent())
fig.text(0.02, 0.27, "Under each panel:\nnumber and % of\nresponsive genes,\nfrom the classification\n(not a label)",
         fontsize=8.5, color="#6B7280", va="center")
fig.savefig("docs/latex/figures/custom_annotated.pdf", bbox_inches="tight")
fig.savefig("docs/figures/custom_annotated.png", dpi=110, bbox_inches="tight")
print("ok")
