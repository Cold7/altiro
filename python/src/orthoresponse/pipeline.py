"""One-call workflow: counts -> five contrasts -> expression models."""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Union

import pandas as pd

from .classify import classify_responses, primed_asymmetry, summarize_classification
from .contrasts import ContrastSet, run_contrasts, write_contrasts

__all__ = ["OrthoResponseResult", "format_asymmetry", "run_orthoresponse", "save_plots"]


@dataclass
class OrthoResponseResult:
    """Result of :func:`run_orthoresponse`."""

    contrasts: ContrastSet
    classification: pd.DataFrame
    summary: pd.DataFrame
    asymmetry: pd.DataFrame

    def __str__(self) -> str:
        """Same text as ``print(res)`` in the R package: summary table and asymmetry test."""
        a, b = self.contrasts.group_names
        summ = self.summary.copy()
        for col in ("pct_of_responsive", "pct_of_all"):
            summ[col] = summ[col].round(1)
        lines = [f"orthoresponse result: A = {a} | B = {b}",
                 f"{len(self.classification)} features classified", "",
                 summ.to_string(index=False, na_rep="NA"), "",
                 format_asymmetry(self.asymmetry)]
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.__str__()

    def save(
        self,
        outdir: Union[str, Path],
        plots: bool = True,
        condition_labels: Sequence[str] = ("C", "S"),
        model_labels: Optional[dict] = None,
        group_names: Optional[Sequence[str]] = None,
        title: Optional[str] = None,
        formats: Sequence[str] = ("png", "pdf"),
        dpi: int = 150,
    ) -> Path:
        """Write all results to ``outdir``.

        Tables (TSV): ``classification``, ``summary``, ``primed_asymmetry``,
        ``normalized_counts`` and ``contrasts/contrast_<key>``. With
        ``plots=True`` (default) also ``model_patterns.<fmt>`` (the 15 patterns
        with counts) and ``model_summary.<fmt>`` (bar chart) for each format in
        ``formats``. ``condition_labels``, ``model_labels``, ``group_names`` and
        ``title`` are passed to :func:`~orthoresponse.plotting.plot_model_patterns`.
        Plots are skipped with a warning if matplotlib is not installed.
        """
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.classification.to_csv(outdir / "classification.tsv", sep="\t")
        self.summary.to_csv(outdir / "summary.tsv", sep="\t", index=False)
        self.asymmetry.to_csv(outdir / "primed_asymmetry.tsv", sep="\t", index=False)
        write_contrasts(self.contrasts, outdir / "contrasts")
        if self.contrasts.normalized_counts is not None:
            self.contrasts.normalized_counts.to_csv(outdir / "normalized_counts.tsv", sep="\t", index_label="feature")
        if plots:
            save_plots(self.classification, outdir, condition_labels=condition_labels, model_labels=model_labels,
                       group_names=group_names, title=title, formats=formats, dpi=dpi)
        return outdir


def format_asymmetry(asymmetry: pd.DataFrame) -> str:
    """One-line description of :func:`~orthoresponse.primed_asymmetry` (same wording as the R package)."""
    a = asymmetry
    return ("Single-group responses ending primed: "
            f"{a['responding_group'].iloc[0]}-only {100 * a['fraction_primed'].iloc[0]:.1f}%, "
            f"{a['responding_group'].iloc[1]}-only {100 * a['fraction_primed'].iloc[1]:.1f}% "
            f"(Fisher p = {a['fisher_pvalue'].iloc[0]:.3g})")


def save_plots(
    classification: pd.DataFrame,
    outdir: Union[str, Path],
    condition_labels: Sequence[str] = ("C", "S"),
    model_labels: Optional[dict] = None,
    group_names: Optional[Sequence[str]] = None,
    title: Optional[str] = None,
    formats: Sequence[str] = ("png", "pdf"),
    dpi: int = 150,
) -> list:
    """Save ``model_patterns.<fmt>`` and ``model_summary.<fmt>`` for a classification.

    Figures are written to files and closed (nothing is shown on screen).
    Returns the written paths; returns an empty list with a warning if
    matplotlib is missing.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        warnings.warn("matplotlib is not installed; plots were not saved. "
                      "Install it with: pip install 'orthoresponse[plot]'", stacklevel=2)
        return []
    from .plotting import plot_model_patterns, plot_summary

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    figs = {
        "model_patterns": plot_model_patterns(classification, group_names=group_names,
                                              condition_labels=condition_labels,
                                              model_labels=model_labels, title=title),
        "model_summary": plot_summary(classification),
    }
    for name, fig in figs.items():
        for fmt in formats:
            path = outdir / f"{name}.{fmt}"
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            written.append(path)
        plt.close(fig)
    return written


def run_orthoresponse(
    counts: pd.DataFrame,
    samples: pd.DataFrame,
    *,
    group_col: str,
    condition_col: str,
    group_a: str,
    group_b: str,
    control: str,
    stress: str,
    preset: str = "paper",
    min_total_count: float = 10,
    alpha: float = 0.05,
    filter_alpha: Optional[float] = None,
    round_counts: bool = False,
    n_cpus: Optional[int] = None,
    **classification_overrides,
) -> OrthoResponseResult:
    """Run PyDESeq2 and classify every feature in one call.

    All arguments are forwarded to :func:`~orthoresponse.run_contrasts` and
    :func:`~orthoresponse.classify_responses`. ``alpha`` is used both for the
    independent filtering (unless ``filter_alpha`` is set) and as the
    classification threshold.

    Examples
    --------
    >>> res = run_orthoresponse(counts, samples, group_col="line", condition_col="treatment",
    ...                         group_a="PL12", group_b="IL22", control="Control", stress="Salt")
    >>> res.summary
    """
    cs = run_contrasts(
        counts, samples, group_col=group_col, condition_col=condition_col, group_a=group_a,
        group_b=group_b, control=control, stress=stress, min_total_count=min_total_count,
        alpha=alpha, filter_alpha=filter_alpha, round_counts=round_counts, n_cpus=n_cpus,
    )
    cls = classify_responses(cs, preset=preset, group_names=(str(group_a), str(group_b)), alpha=alpha,
                             **classification_overrides)
    return OrthoResponseResult(contrasts=cs, classification=cls, summary=summarize_classification(cls),
                               asymmetry=primed_asymmetry(cls))
