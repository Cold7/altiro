"""Differential-expression contrasts for the 2 x 2 design.

:func:`run_contrasts` fits one negative-binomial model with PyDESeq2 using a
four-level factor (``A_control``, ``A_stress``, ``B_control``, ``B_stress``)
and extracts the five contrasts required by the classifier. This mirrors the
DESeq2 script of the original publication (single ``~ condition`` design over
all samples, size factors estimated on both groups together).

:func:`read_contrasts` loads contrast tables produced elsewhere (e.g. by
DESeq2 in R, edgeR or limma exported to the same column names).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Mapping, Optional, Union

import numpy as np
import pandas as pd

from .classify import CONTRAST_KEYS

__all__ = ["ContrastSet", "build_design", "run_contrasts", "read_contrasts", "write_contrasts"]

CELLS = ("A_control", "A_stress", "B_control", "B_stress")

_CONTRAST_CELLS = {
    "a_response": ("A_stress", "A_control"),
    "b_response": ("B_stress", "B_control"),
    "basal": ("B_control", "A_control"),
    "b_stress_vs_a_control": ("B_stress", "A_control"),
    "a_stress_vs_b_control": ("A_stress", "B_control"),
}


@dataclass
class ContrastSet:
    """Container for the five contrast tables and the normalized counts."""

    tables: Dict[str, pd.DataFrame]
    normalized_counts: Optional[pd.DataFrame] = None
    design: Optional[pd.DataFrame] = None
    group_names: tuple = ("A", "B")
    info: dict = field(default_factory=dict)

    def __getitem__(self, key: str) -> pd.DataFrame:
        return self.tables[key]

    def keys(self):
        return self.tables.keys()


def build_design(
    samples: pd.DataFrame,
    group_col: str,
    condition_col: str,
    group_a: str,
    group_b: str,
    control: str,
    stress: str,
) -> pd.DataFrame:
    """Map sample annotations to the four design cells.

    Samples whose group or condition is not one of the four requested levels
    are dropped with a warning. Raises if any cell is empty.
    """
    for col in (group_col, condition_col):
        if col not in samples.columns:
            raise KeyError(f"Column {col!r} not found in the sample table")
    g = samples[group_col].astype(str)
    c = samples[condition_col].astype(str)
    gmap = {str(group_a): "A", str(group_b): "B"}
    cmap = {str(control): "control", str(stress): "stress"}
    keep = g.isin(gmap) & c.isin(cmap)
    if (~keep).any():
        warnings.warn(f"Dropping {int((~keep).sum())} sample(s) outside the 2 x 2 design: "
                      f"{list(samples.index[~keep])}", stacklevel=2)
    d = pd.DataFrame(
        {"group": g[keep], "condition": c[keep]},
        index=samples.index[keep],
    )
    d["cell"] = g[keep].map(gmap) + "_" + c[keep].map(cmap)
    sizes = d["cell"].value_counts().reindex(CELLS, fill_value=0)
    if (sizes == 0).any():
        raise ValueError(f"Empty design cell(s): {list(sizes.index[sizes == 0])}. Sizes: {sizes.to_dict()}")
    if (sizes < 2).any():
        warnings.warn(f"Some design cells have fewer than 2 replicates: {sizes.to_dict()}", stacklevel=2)
    return d


def run_contrasts(
    counts: pd.DataFrame,
    samples: pd.DataFrame,
    *,
    group_col: str,
    condition_col: str,
    group_a: str,
    group_b: str,
    control: str,
    stress: str,
    min_total_count: float = 10,
    alpha: float = 0.05,
    filter_alpha: Optional[float] = None,
    cooks_filter: bool = True,
    independent_filter: bool = True,
    round_counts: bool = False,
    n_cpus: Optional[int] = None,
    quiet: bool = True,
) -> ContrastSet:
    """Fit PyDESeq2 once and return the five contrasts.

    Parameters
    ----------
    counts:
        Raw counts, features in rows and samples in columns.
    samples:
        Sample annotation indexed by sample name (matching ``counts`` columns).
    group_col, condition_col:
        Columns of ``samples`` holding the group and the condition.
    group_a, group_b:
        Levels of ``group_col``. A is the reference (e.g. the sensitive
        species), B the comparison (e.g. the tolerant species). In the original
        publication A = *S. lycopersicum* and B = *S. pennellii*.
    control, stress:
        Levels of ``condition_col``.
    min_total_count:
        Features with a total count ``<=`` this value are removed before fitting
        (the published script used ``rowSums > 10``).
    alpha:
        Significance level passed to the independent filtering unless
        ``filter_alpha`` is given. The published script used 0.01 for the
        filtering and 0.05 for the classification.
    round_counts:
        Round non-integer input (e.g. estimated counts) instead of failing.

    Returns
    -------
    ContrastSet
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as err:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "run_contrasts() needs PyDESeq2. Install it with: pip install 'orthoresponse[deseq]'"
        ) from err

    design = build_design(samples, group_col, condition_col, group_a, group_b, control, stress)
    missing = [s for s in design.index if s not in counts.columns]
    if missing:
        raise KeyError(f"Samples missing from the count matrix: {missing}")
    mat = counts.loc[:, design.index].apply(pd.to_numeric, errors="raise")
    if mat.isna().any().any():
        raise ValueError("The count matrix contains missing values")
    if not np.allclose(mat.to_numpy(), np.round(mat.to_numpy())):
        if not round_counts:
            raise ValueError("Counts are not integers. Use raw counts, or set round_counts=True for estimated counts.")
        mat = mat.round()
    if (mat.to_numpy() < 0).any():
        raise ValueError("Counts must be non-negative")
    mat = mat.astype(np.int64)
    keep = mat.sum(axis=1) > min_total_count
    mat = mat.loc[keep]
    if mat.index.has_duplicates:
        raise ValueError("Duplicated feature identifiers in the count matrix")

    meta = design[["cell"]].copy()
    dds = DeseqDataSet(counts=mat.T, metadata=meta, design="~cell", quiet=quiet, n_cpus=n_cpus)
    dds.deseq2()

    f_alpha = alpha if filter_alpha is None else filter_alpha
    tables = {}
    for key in CONTRAST_KEYS:
        num, den = _CONTRAST_CELLS[key]
        st = DeseqStats(dds, contrast=["cell", num, den], alpha=f_alpha, cooks_filter=cooks_filter,
                        independent_filter=independent_filter, quiet=quiet, n_cpus=n_cpus)
        st.summary()
        res = st.results_df.copy()
        res.index.name = "feature"
        tables[key] = res

    sf = np.asarray(dds.obs["size_factors"]) if "size_factors" in dds.obs else np.asarray(dds.obsm["size_factors"])
    norm_counts = pd.DataFrame(mat.to_numpy() / sf[np.newaxis, :], index=mat.index, columns=mat.columns)
    info = {
        "n_input_features": int(len(counts)),
        "n_tested_features": int(len(mat)),
        "min_total_count": min_total_count,
        "filter_alpha": f_alpha,
        "engine": "PyDESeq2",
    }
    return ContrastSet(tables=tables, normalized_counts=norm_counts, design=design,
                       group_names=(str(group_a), str(group_b)), info=info)


def read_contrasts(
    paths: Mapping[str, Union[str, Path]],
    id_col: Optional[str] = None,
    sep: Optional[str] = None,
    lfc_col: str = "log2FoldChange",
    padj_col: str = "padj",
    se_col: str = "lfcSE",
) -> Dict[str, pd.DataFrame]:
    """Read five contrast tables from delimited files.

    Parameters
    ----------
    paths:
        Mapping from contrast key (see :data:`~orthoresponse.CONTRAST_KEYS`) to a
        file path. For the files written by the original publication script use::

            {"a_response": "Sly.drought_Sly.control.csv",
             "b_response": "Spenn.drought_Spenn.control.csv",
             "basal": "Spenn.control_Sly.control.csv",
             "b_stress_vs_a_control": "Spenn.drought_Sly.control.csv",
             "a_stress_vs_b_control": "Sly.drought_Spenn.control.csv"}

        with ``id_col="Row.names"``.
    id_col:
        Column holding the feature identifier. Defaults to the first column.
    sep:
        Field separator; inferred from the extension when omitted.
    lfc_col, padj_col, se_col:
        Column names to map onto ``log2FoldChange``, ``padj`` and ``lfcSE``.
    """
    out = {}
    for key in CONTRAST_KEYS:
        if key not in paths:
            raise KeyError(f"No path given for contrast {key!r}")
        p = Path(paths[key])
        s = sep or ("\t" if p.suffix.lower() in (".tsv", ".txt", ".tab") or p.name.endswith((".tsv.gz", ".txt.gz")) else ",")
        df = pd.read_csv(p, sep=s)
        col = id_col or df.columns[0]
        if col not in df.columns:
            raise KeyError(f"{p}: id column {col!r} not found")
        df = df.set_index(col)
        df.index = df.index.astype(str)
        df.index.name = "feature"
        df = df.rename(columns={lfc_col: "log2FoldChange", padj_col: "padj", se_col: "lfcSE"})
        out[key] = df
    return out


def write_contrasts(contrasts, outdir: Union[str, Path], sep: str = "\t") -> Dict[str, Path]:
    """Write each contrast table to ``<outdir>/contrast_<key>.tsv``."""
    tables = contrasts.tables if hasattr(contrasts, "tables") else contrasts
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key in CONTRAST_KEYS:
        p = outdir / f"contrast_{key}.tsv"
        tables[key].to_csv(p, sep=sep, index_label="feature")
        paths[key] = p
    return paths
