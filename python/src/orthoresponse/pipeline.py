"""One-call workflow: counts -> five contrasts -> expression models."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from .classify import classify_responses, primed_asymmetry, summarize_classification
from .contrasts import ContrastSet, run_contrasts, write_contrasts

__all__ = ["OrthoResponseResult", "run_orthoresponse"]


@dataclass
class OrthoResponseResult:
    """Result of :func:`run_orthoresponse`."""

    contrasts: ContrastSet
    classification: pd.DataFrame
    summary: pd.DataFrame
    asymmetry: pd.DataFrame

    def save(self, outdir: Union[str, Path]) -> Path:
        """Write classification, summary, asymmetry test and contrast tables as TSV files."""
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.classification.to_csv(outdir / "classification.tsv", sep="\t")
        self.summary.to_csv(outdir / "summary.tsv", sep="\t", index=False)
        self.asymmetry.to_csv(outdir / "primed_asymmetry.tsv", sep="\t", index=False)
        write_contrasts(self.contrasts, outdir / "contrasts")
        if self.contrasts.normalized_counts is not None:
            self.contrasts.normalized_counts.to_csv(outdir / "normalized_counts.tsv", sep="\t", index_label="feature")
        return outdir


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
