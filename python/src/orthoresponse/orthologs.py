"""Ortholog helpers: reciprocal best BLAST hits and merged count matrices.

These functions cover the cross-species part of the workflow used in the
original publication: orthologs defined as reciprocal best BLASTp hits, and a
single count matrix in which each row is an ortholog pair and the columns are
the samples of both species.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

__all__ = [
    "BLAST6_COLUMNS",
    "read_blast_tabular",
    "best_hits",
    "reciprocal_best_hits",
    "merge_ortholog_counts",
]

BLAST6_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]

PathOrFrame = Union[str, Path, pd.DataFrame]


def read_blast_tabular(path: Union[str, Path], strip_pattern: Optional[str] = None) -> pd.DataFrame:
    """Read BLAST ``-outfmt 6`` output (12 standard columns).

    Parameters
    ----------
    strip_pattern:
        Optional regular expression removed from query and subject identifiers,
        e.g. ``r"\\.\\d+$"`` to turn ``Solyc01g005000.3.1`` into
        ``Solyc01g005000.3`` or ``r"\\.\\d+\\.\\d+$"`` to get the gene id.
    """
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 12:
        raise ValueError(f"{path}: expected at least 12 columns (BLAST -outfmt 6), found {df.shape[1]}")
    df = df.iloc[:, :12]
    df.columns = BLAST6_COLUMNS
    df["qseqid"] = df["qseqid"].astype(str)
    df["sseqid"] = df["sseqid"].astype(str)
    if strip_pattern:
        rx = re.compile(strip_pattern)
        df["qseqid"] = df["qseqid"].map(lambda x: rx.sub("", x))
        df["sseqid"] = df["sseqid"].map(lambda x: rx.sub("", x))
    return df


def _load(x: PathOrFrame, strip_pattern: Optional[str]) -> pd.DataFrame:
    if isinstance(x, pd.DataFrame):
        df = x.copy()
        if strip_pattern:
            rx = re.compile(strip_pattern)
            df["qseqid"] = df["qseqid"].astype(str).map(lambda v: rx.sub("", v))
            df["sseqid"] = df["sseqid"].astype(str).map(lambda v: rx.sub("", v))
        return df
    return read_blast_tabular(x, strip_pattern)


def best_hits(blast: pd.DataFrame) -> pd.DataFrame:
    """Best subject per query: highest bitscore, then lowest e-value, then subject id.

    Self-identical query/subject pairs are kept (they are meaningful across
    species). Ties are broken deterministically by subject identifier.
    """
    df = blast.sort_values(["qseqid", "bitscore", "evalue", "sseqid"], ascending=[True, False, True, True])
    return df.drop_duplicates("qseqid", keep="first").reset_index(drop=True)


def reciprocal_best_hits(
    blast_ab: PathOrFrame,
    blast_ba: PathOrFrame,
    *,
    max_evalue: Optional[float] = None,
    min_bitscore: Optional[float] = None,
    strip_pattern: Optional[str] = None,
) -> pd.DataFrame:
    """Orthologs as reciprocal best hits between two proteomes.

    Parameters
    ----------
    blast_ab:
        BLAST of proteome A (queries) against proteome B (database).
    blast_ba:
        BLAST of proteome B against proteome A.
    max_evalue, min_bitscore:
        Optional filters applied before choosing best hits.
    strip_pattern:
        Regular expression removed from identifiers (see :func:`read_blast_tabular`).

    Returns
    -------
    DataFrame with columns ``ortholog_id``, ``gene_a``, ``gene_b``,
    ``pident_ab``, ``evalue_ab``, ``bitscore_ab``, ``evalue_ba``, ``bitscore_ba``.
    If several isoforms of the same gene are collapsed by ``strip_pattern``, the
    best-scoring isoform pair is kept.
    """
    ab = _load(blast_ab, strip_pattern)
    ba = _load(blast_ba, strip_pattern)
    for df in (ab, ba):
        if max_evalue is not None:
            df.drop(df.index[df["evalue"] > max_evalue], inplace=True)
        if min_bitscore is not None:
            df.drop(df.index[df["bitscore"] < min_bitscore], inplace=True)
    bab = best_hits(ab)[["qseqid", "sseqid", "pident", "evalue", "bitscore"]]
    bba = best_hits(ba)[["qseqid", "sseqid", "evalue", "bitscore"]]
    m = bab.merge(bba, left_on=["qseqid", "sseqid"], right_on=["sseqid", "qseqid"], suffixes=("_ab", "_ba"))
    out = pd.DataFrame(
        {
            "gene_a": m["qseqid_ab"],
            "gene_b": m["sseqid_ab"],
            "pident_ab": m["pident"],
            "evalue_ab": m["evalue_ab"],
            "bitscore_ab": m["bitscore_ab"],
            "evalue_ba": m["evalue_ba"],
            "bitscore_ba": m["bitscore_ba"],
        }
    )
    out = out.sort_values("gene_a").reset_index(drop=True)
    out.insert(0, "ortholog_id", [f"ortholog_{i + 1}" for i in range(len(out))])
    return out


def merge_ortholog_counts(
    counts_a: pd.DataFrame,
    counts_b: pd.DataFrame,
    orthologs: pd.DataFrame,
    *,
    gene_a_col: str = "gene_a",
    gene_b_col: str = "gene_b",
    id_col: Optional[str] = "ortholog_id",
    lengths_a: Optional[pd.Series] = None,
    lengths_b: Optional[pd.Series] = None,
    length_normalize: bool = False,
) -> pd.DataFrame:
    """Combine two species' count matrices into one ortholog-level matrix.

    Parameters
    ----------
    counts_a, counts_b:
        Raw counts (features x samples) for species A and B. Sample names must
        be unique across both matrices.
    orthologs:
        One-to-one ortholog table (e.g. from :func:`reciprocal_best_hits`).
    id_col:
        Column of ``orthologs`` used as row name. If None, rows are named
        ``"<gene_a>|<gene_b>"``.
    lengths_a, lengths_b:
        Transcript/gene lengths in base pairs, indexed by gene id. Required when
        ``length_normalize`` is True.
    length_normalize:
        Apply the length correction used in the original publication,
        ``Q = round(E / L_kb)``, to each species before merging. This puts
        orthologs of different lengths on a comparable scale but changes the
        count distribution assumed by DESeq2; keep it False unless you want to
        replicate that analysis.

    Returns
    -------
    DataFrame with one row per ortholog pair and the samples of A followed by B.
    """
    overlap = set(counts_a.columns) & set(counts_b.columns)
    if overlap:
        raise ValueError(f"Sample names shared by both matrices: {sorted(overlap)[:5]}")
    orth = orthologs[[c for c in (id_col, gene_a_col, gene_b_col) if c]].copy()
    if orth[gene_a_col].duplicated().any() or orth[gene_b_col].duplicated().any():
        raise ValueError("The ortholog table must be one-to-one (duplicated gene ids found)")
    present = orth[gene_a_col].isin(counts_a.index) & orth[gene_b_col].isin(counts_b.index)
    if (~present).any():
        import warnings

        warnings.warn(f"{int((~present).sum())} ortholog pair(s) absent from the count matrices were dropped",
                      stacklevel=2)
    orth = orth[present]
    a = counts_a.loc[orth[gene_a_col]].astype(float)
    b = counts_b.loc[orth[gene_b_col]].astype(float)
    if length_normalize:
        if lengths_a is None or lengths_b is None:
            raise ValueError("length_normalize=True requires lengths_a and lengths_b")
        la = pd.Series(lengths_a).reindex(orth[gene_a_col]).to_numpy(dtype=float) / 1000.0
        lb = pd.Series(lengths_b).reindex(orth[gene_b_col]).to_numpy(dtype=float) / 1000.0
        if np.isnan(la).any() or np.isnan(lb).any():
            raise ValueError("Missing lengths for some orthologs")
        a = np.round(a.div(la, axis=0))
        b = np.round(b.div(lb, axis=0))
    rows = orth[id_col].astype(str).to_numpy() if id_col else (orth[gene_a_col] + "|" + orth[gene_b_col]).to_numpy()
    merged = pd.concat([a.reset_index(drop=True), b.reset_index(drop=True)], axis=1)
    merged.index = rows
    merged.index.name = "feature"
    if np.allclose(merged.to_numpy(), np.round(merged.to_numpy())):
        merged = merged.round().astype(np.int64)
    return merged
