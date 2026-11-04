"""Plots for expression models (requires matplotlib: ``pip install 'orthoresponse[plot]'``)."""

from __future__ import annotations

from typing import Optional, Sequence

import pandas as pd

from .classify import summarize_classification

__all__ = ["PATTERNS", "plot_model_patterns", "plot_summary", "plot_feature"]

COLOR_A = "#1f77b4"
COLOR_B = "#ff7f0e"

# Schematic levels (A control, A stress, B control, B stress) for each model/case.
PATTERNS = {
    ("stress_ready", 1): (0, 0, -1, 0),
    ("stress_ready", 2): (0, 0, 1, 0),
    ("stress_ready", 3): (0, 1, 1, 1),
    ("stress_ready", 4): (0, -1, -1, -1),
    ("shared", 1): (0, 1, 0, 1),
    ("shared", 2): (0, -1, 0, -1),
    ("unique", 1): (0, 1, 0, 0),
    ("unique", 2): (0, -1, 0, 0),
    ("unique", 3): (0, 0, 0, 1),
    ("unique", 4): (0, 0, 0, -1),
    ("opposite", 1): (0, 1, 0, -1),
    ("opposite", 2): (0, -1, 0, 1),
    ("no_response", 1): (0, 0, 1, 1),
    ("no_response", 2): (0, 0, -1, -1),
    ("no_response", 3): (0, 0, 0, 0),
}

_ROWS = [
    ("stress_ready", [1, 2, 3, 4]),
    ("shared", [1, 2]),
    ("unique", [1, 2, 3, 4]),
    ("opposite", [1, 2]),
    ("no_response", [1, 2, 3]),
]

_TITLES = {
    "stress_ready": "Stress-ready",
    "shared": "Shared response",
    "unique": "Unique response",
    "opposite": "Opposite response",
    "no_response": "No response",
}


def _mini(ax, levels):
    a_c, a_s, b_c, b_s = levels
    ax.plot([0, 1, 2, 3], [a_c, a_c, a_s, a_s], color=COLOR_A, lw=2.5)
    ax.plot([3, 4, 5, 6], [b_c, b_c, b_s, b_s], color=COLOR_B, lw=2.5)
    for x, ls in ((1.5, "--"), (3, "-"), (4.5, "--")):
        ax.axvline(x, color="k", ls=ls, alpha=0.25, lw=1)
    ax.axhline(0, color="k", ls="--", alpha=0.25, lw=1)
    ax.set_ylim(-1.3, 1.3)
    ax.set_xlim(0, 6)
    ax.set_xticks([0.75, 2.25, 3.75, 5.25])
    ax.set_xticklabels(["C", "S", "C", "S"], fontsize=7)
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, length=0)
    ax.set_yticks([])


def plot_model_patterns(
    classification: Optional[pd.DataFrame] = None,
    group_names: Optional[Sequence[str]] = None,
    figsize=(11, 10),
):
    """Draw the 15 expression patterns, annotated with counts if a classification is given.

    C = control, S = stress. The first half of each panel is group A, the
    second half group B.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    if group_names is None:
        group_names = classification.attrs.get("group_names", ("A", "B")) if classification is not None else ("A", "B")
    summ = summarize_classification(classification).set_index(["model", "case"]) if classification is not None else None

    fig, axes = plt.subplots(len(_ROWS), 4, figsize=figsize)
    for r, (model, cases) in enumerate(_ROWS):
        for c in range(4):
            ax = axes[r, c]
            if c >= len(cases):
                ax.axis("off")
                continue
            case = cases[c]
            _mini(ax, PATTERNS[(model, case)])
            label = f"Case {case}"
            if summ is not None:
                n = int(summ["n"].get((model, case), 0))
                pct = summ["pct_of_responsive"].get((model, case), float("nan"))
                label += f": {n}" + (f" ({pct:.1f}%)" if pd.notna(pct) else "")
            ax.set_xlabel(label, fontsize=8)
        axes[r, 0].set_ylabel(_TITLES[model], fontsize=9)
    handles = [Line2D([0], [0], color=COLOR_A, lw=3), Line2D([0], [0], color=COLOR_B, lw=3)]
    fig.legend(handles, list(group_names), loc="lower right", bbox_to_anchor=(0.97, 0.08), frameon=False)
    fig.suptitle("Expression models (C = control, S = stress); % of responsive features", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


def plot_summary(classification: pd.DataFrame, ax=None):
    """Horizontal bar chart of the number of features per model and case."""
    import matplotlib.pyplot as plt

    s = summarize_classification(classification)
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 0.3 * len(s) + 1))
    ax.barh(s["label"][::-1], s["n"][::-1], color="#4c72b0")
    for y, v in enumerate(s["n"][::-1]):
        ax.text(v, y, f" {v}", va="center", fontsize=8)
    ax.set_xlabel("Features")
    ax.figure.tight_layout()
    return ax.figure


def plot_feature(normalized_counts: pd.DataFrame, design: pd.DataFrame, feature: str, ax=None, log=True):
    """Dot plot of one feature's normalized counts across the four design cells."""
    import matplotlib.pyplot as plt
    import numpy as np

    if ax is None:
        _, ax = plt.subplots(figsize=(4, 3))
    vals = normalized_counts.loc[feature, design.index]
    if log:
        vals = np.log2(vals + 1)
    cells = ["A_control", "A_stress", "B_control", "B_stress"]
    for i, cell in enumerate(cells):
        v = vals[design["cell"] == cell]
        color = COLOR_A if cell.startswith("A") else COLOR_B
        ax.scatter(np.full(len(v), i), v, color=color, zorder=3)
        ax.hlines(v.mean(), i - 0.25, i + 0.25, color="k")
    ax.set_xticks(range(4))
    ax.set_xticklabels(cells, rotation=30)
    ax.set_ylabel("log2(normalized count + 1)" if log else "normalized count")
    ax.set_title(str(feature))
    ax.figure.tight_layout()
    return ax.figure
