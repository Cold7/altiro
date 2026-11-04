"""Rule-based classification of expression response models in a 2 x 2 design.

The design has two *groups* (species, genotypes, cell lines, ...) called A and B,
each measured under a *control* and a *stress* condition. Every feature (an
ortholog pair, or a gene when both groups share a genome) is assigned to one
expression model by a decision tree that uses five differential-expression
contrasts:

=========================  ============================================
key                        contrast (numerator vs denominator)
=========================  ============================================
``a_response``             A stress  vs A control
``b_response``             B stress  vs B control
``basal``                  B control vs A control
``b_stress_vs_a_control``  B stress  vs A control
``a_stress_vs_b_control``  A stress  vs B control
=========================  ============================================

With the default ``preset="paper"`` the rules reproduce exactly the decision
tree published in Contreras-Riquelme et al. (2025, PLoS One 20(5): e0324724),
where A = *Solanum lycopersicum* and B = *Solanum pennellii*.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm

__all__ = [
    "CONTRAST_KEYS",
    "CONTRAST_LABELS",
    "ClassificationOptions",
    "classification_options",
    "model_definitions",
    "equivalence_test",
    "adjust_bh",
    "classify_responses",
    "summarize_classification",
    "primed_asymmetry",
]

CONTRAST_KEYS: Tuple[str, ...] = (
    "a_response",
    "b_response",
    "basal",
    "b_stress_vs_a_control",
    "a_stress_vs_b_control",
)

CONTRAST_LABELS: Dict[str, str] = {
    "a_response": "A stress vs A control",
    "b_response": "B stress vs B control",
    "basal": "B control vs A control",
    "b_stress_vs_a_control": "B stress vs A control",
    "a_stress_vs_b_control": "A stress vs B control",
}

RESPONSIVE_MODELS = ("stress_ready", "shared", "unique", "opposite", "unresolved")

_NA_ACTIONS = ("not_de", "unclassified")


@dataclass(frozen=True)
class ClassificationOptions:
    """Options controlling how differential-expression calls are turned into models.

    Attributes
    ----------
    alpha:
        Adjusted p-value threshold for a significant difference (default 0.05).
    min_abs_lfc:
        Minimum absolute log2 fold change required, in addition to ``alpha``,
        to call a contrast significant (default 0, i.e. sign only).
    na_action:
        ``"not_de"`` (paper behaviour): a missing adjusted p-value (independent
        filtering, Cook's outliers) counts as *not differentially expressed*.
        ``"unclassified"``: features with a missing value in any contrast that
        the decision tree consults are labelled ``unclassified``.
    require_basal_difference:
        If True, a stress-ready call also requires a significant basal
        difference (``basal`` contrast) in the direction implied by the case.
        Features failing this check are labelled as the corresponding unique
        response.
    equivalence_margin:
        If set (log2 units, e.g. 0.5), "same expression level" must be
        demonstrated with a two one-sided test (TOST) on the matching contrast
        instead of being inferred from a non-significant difference. Features
        that are neither different nor equivalent are labelled ``unresolved``.
        Requires an ``lfcSE`` column.
    sign_digits:
        Number of decimals used to round log2 fold changes before taking their
        sign (the published script used 2).
    """

    alpha: float = 0.05
    min_abs_lfc: float = 0.0
    na_action: str = "not_de"
    require_basal_difference: bool = False
    equivalence_margin: Optional[float] = None
    sign_digits: int = 2

    def __post_init__(self) -> None:
        if not 0 < self.alpha < 1:
            raise ValueError("alpha must be in (0, 1)")
        if self.min_abs_lfc < 0:
            raise ValueError("min_abs_lfc must be >= 0")
        if self.na_action not in _NA_ACTIONS:
            raise ValueError(f"na_action must be one of {_NA_ACTIONS}")
        if self.equivalence_margin is not None and self.equivalence_margin <= 0:
            raise ValueError("equivalence_margin must be > 0 or None")


_PRESETS: Dict[str, ClassificationOptions] = {
    "paper": ClassificationOptions(),
    "strict": ClassificationOptions(na_action="unclassified", require_basal_difference=True),
}


def classification_options(preset: str = "paper", **overrides) -> ClassificationOptions:
    """Return the options of a preset, optionally overriding individual fields.

    Presets
    -------
    ``"paper"``
        Exact rules of the published decision tree.
    ``"strict"``
        ``na_action="unclassified"`` and ``require_basal_difference=True``.

    Examples
    --------
    >>> classification_options("paper", equivalence_margin=0.5).equivalence_margin
    0.5
    """
    if preset not in _PRESETS:
        raise ValueError(f"Unknown preset {preset!r}; choose from {sorted(_PRESETS)}")
    unknown = set(overrides) - set(ClassificationOptions.__dataclass_fields__)
    if unknown:
        raise TypeError(f"Unknown option(s): {sorted(unknown)}")
    overrides = {k: v for k, v in overrides.items() if v is not None}
    return replace(_PRESETS[preset], **overrides)


# ----------------------------------------------------------------------------
# Model definitions
# ----------------------------------------------------------------------------

_DEFINITIONS = [
    # model, case, primed, responding, description template
    ("stress_ready", 1, "A", "B", "{B} is induced under stress up to the control level of {A} ({A} is primed)"),
    ("stress_ready", 2, "A", "B", "{B} is repressed under stress down to the control level of {A} ({A} is primed)"),
    ("stress_ready", 3, "B", "A", "{A} is induced under stress up to the control level of {B} ({B} is primed)"),
    ("stress_ready", 4, "B", "A", "{A} is repressed under stress down to the control level of {B} ({B} is primed)"),
    ("shared", 1, "", "both", "Induced under stress in both {A} and {B}"),
    ("shared", 2, "", "both", "Repressed under stress in both {A} and {B}"),
    ("unique", 1, "", "A", "Induced under stress only in {A}"),
    ("unique", 2, "", "A", "Repressed under stress only in {A}"),
    ("unique", 3, "", "B", "Induced under stress only in {B}"),
    ("unique", 4, "", "B", "Repressed under stress only in {B}"),
    ("opposite", 1, "", "both", "Induced in {A} and repressed in {B}"),
    ("opposite", 2, "", "both", "Repressed in {A} and induced in {B}"),
    ("no_response", 1, "", "none", "No response; basal level higher in {B} than in {A}"),
    ("no_response", 2, "", "none", "No response; basal level lower in {B} than in {A}"),
    ("no_response", 3, "", "none", "No response; no basal difference detected"),
    ("unresolved", 1, "", "B", "Responds only in {B}; level under stress neither different from nor equivalent to {A} control"),
    ("unresolved", 2, "", "B", "Responds only in {B}; level under stress neither different from nor equivalent to {A} control"),
    ("unresolved", 3, "", "A", "Responds only in {A}; level under stress neither different from nor equivalent to {B} control"),
    ("unresolved", 4, "", "A", "Responds only in {A}; level under stress neither different from nor equivalent to {B} control"),
    ("unclassified", 0, "", "", "Missing statistics in a contrast required by the decision tree"),
]


def model_definitions(group_names: Sequence[str] = ("A", "B"), include_extra: bool = False) -> pd.DataFrame:
    """Table describing every expression model and case.

    Parameters
    ----------
    group_names:
        Display names for groups A and B, used in the descriptions.
    include_extra:
        Also list ``unresolved`` and ``unclassified``, which are only produced
        by non-default options.
    """
    a, b = group_names
    rows = []
    for model, case, primed, responding, text in _DEFINITIONS:
        if not include_extra and model in ("unresolved", "unclassified"):
            continue
        rows.append(
            {
                "model": model,
                "case": case,
                "label": model if model == "unclassified" else f"{model}_{case}",
                "primed_group": {"A": a, "B": b}.get(primed, ""),
                "responding_group": {"A": a, "B": b}.get(responding, responding),
                "description": text.format(A=a, B=b),
            }
        )
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Statistics helpers
# ----------------------------------------------------------------------------

def adjust_bh(p: Union[Sequence[float], np.ndarray]) -> np.ndarray:
    """Benjamini-Hochberg adjustment that ignores missing values (like R's ``p.adjust``)."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return out
    pv = p[ok]
    order = np.argsort(pv)[::-1]
    ranks = np.arange(m, 0, -1)
    adj = np.minimum.accumulate(pv[order] * m / ranks)
    res = np.empty(m)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


def equivalence_test(log2fc, lfc_se, margin: float) -> pd.DataFrame:
    """Two one-sided Wald tests (TOST) for ``|log2FC| < margin``.

    Returns a DataFrame with ``pvalue`` (the larger of the two one-sided
    p-values) and ``padj`` (Benjamini-Hochberg). This is the same Wald-based
    construction used by DESeq2 with ``altHypothesis = "lessAbs"``.
    """
    lfc = np.asarray(log2fc, dtype=float)
    se = np.asarray(lfc_se, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p_low = norm.sf((lfc + margin) / se)   # H0: lfc <= -margin
        p_high = norm.cdf((lfc - margin) / se)  # H0: lfc >= +margin
    p = np.fmax(p_low, p_high)
    p[np.isnan(lfc) | np.isnan(se) | (se <= 0)] = np.nan
    idx = getattr(log2fc, "index", None)
    return pd.DataFrame({"pvalue": p, "padj": adjust_bh(p)}, index=idx)


# ----------------------------------------------------------------------------
# Input handling
# ----------------------------------------------------------------------------

ContrastInput = Union[Mapping[str, pd.DataFrame], pd.DataFrame, "object"]


def _as_tables(contrasts) -> Dict[str, pd.DataFrame]:
    """Accept a mapping of five DataFrames, a ContrastSet, or one wide DataFrame."""
    if hasattr(contrasts, "tables"):
        contrasts = contrasts.tables
    if isinstance(contrasts, pd.DataFrame):
        tables = {}
        for key in CONTRAST_KEYS:
            cols = {c: c[len(key) + 1:] for c in contrasts.columns if c.startswith(key + "_")}
            if not cols:
                raise KeyError(f"No columns starting with '{key}_' in the wide table")
            tables[key] = contrasts[list(cols)].rename(columns=cols)
        contrasts = tables
    missing = [k for k in CONTRAST_KEYS if k not in contrasts]
    if missing:
        raise KeyError(f"Missing contrasts: {missing}. Required keys: {list(CONTRAST_KEYS)}")
    out = {}
    for key in CONTRAST_KEYS:
        df = contrasts[key]
        for col in ("log2FoldChange", "padj"):
            if col not in df.columns:
                raise KeyError(f"Contrast '{key}' lacks column '{col}'")
        if df.index.has_duplicates:
            raise ValueError(f"Contrast '{key}' has duplicated feature identifiers")
        out[key] = df
    return out


def _align(tables: Dict[str, pd.DataFrame]) -> Tuple[pd.Index, Dict[str, pd.DataFrame]]:
    index = tables[CONTRAST_KEYS[0]].index
    for key in CONTRAST_KEYS[1:]:
        index = index.intersection(tables[key].index, sort=False)
    sizes = {k: len(v) for k, v in tables.items()}
    if len(set(sizes.values())) > 1 or len(index) != sizes[CONTRAST_KEYS[0]]:
        import warnings

        warnings.warn(
            f"Contrast tables cover different features {sizes}; using the {len(index)} shared ones.",
            stacklevel=3,
        )
    return index, {k: v.loc[index] for k, v in tables.items()}


# ----------------------------------------------------------------------------
# Classification
# ----------------------------------------------------------------------------

def _sign(lfc: np.ndarray, digits: int) -> np.ndarray:
    r = np.round(lfc, digits)
    with np.errstate(invalid="ignore"):
        return np.where(r > 0, 1, np.where(r < 0, -1, 0))


def classify_responses(
    contrasts,
    preset: str = "paper",
    group_names: Sequence[str] = ("A", "B"),
    options: Optional[ClassificationOptions] = None,
    keep_statistics: bool = True,
    **overrides,
) -> pd.DataFrame:
    """Assign every feature to an expression model.

    Parameters
    ----------
    contrasts:
        Either a mapping ``{key: DataFrame}`` with the five keys in
        :data:`CONTRAST_KEYS` (each DataFrame indexed by feature id with at
        least ``log2FoldChange`` and ``padj``; ``lfcSE`` is needed for
        equivalence tests), a :class:`~orthoresponse.contrasts.ContrastSet`, or
        one wide DataFrame whose columns are ``<key>_log2FoldChange``,
        ``<key>_padj`` ...
    preset:
        ``"paper"`` (default) or ``"strict"``; see :func:`classification_options`.
    group_names:
        Display names for groups A and B (e.g. ``("S. lycopersicum", "S. pennellii")``).
    options:
        A :class:`ClassificationOptions` instance; overrides ``preset``.
    keep_statistics:
        Append ``<key>_log2FoldChange`` and ``<key>_padj`` columns.
    **overrides:
        Any field of :class:`ClassificationOptions`, e.g. ``alpha=0.01``.

    Returns
    -------
    pandas.DataFrame
        Indexed by feature with columns ``label``, ``model``, ``case``,
        ``primed_group``, ``responding_group``, ``description`` and ``flags``.
    """
    opts = options if options is not None else classification_options(preset, **overrides)
    tables = _as_tables(contrasts)
    index, t = _align(tables)
    n = len(index)

    def stats(key):
        df = t[key]
        padj = pd.to_numeric(df["padj"], errors="coerce").to_numpy(dtype=float)
        lfc = pd.to_numeric(df["log2FoldChange"], errors="coerce").to_numpy(dtype=float)
        na = np.isnan(padj)
        with np.errstate(invalid="ignore"):
            de = (~na) & (padj < opts.alpha) & (np.abs(np.nan_to_num(lfc)) >= opts.min_abs_lfc)
        return de, na, _sign(lfc, opts.sign_digits), lfc

    de_a, na_a, s_a, _ = stats("a_response")
    de_b, na_b, s_b, _ = stats("b_response")
    de_bas, na_bas, s_bas, _ = stats("basal")
    de_mb, na_mb, _, _ = stats("b_stress_vs_a_control")
    de_ma, na_ma, _, _ = stats("a_stress_vs_b_control")

    if opts.equivalence_margin is not None:
        for key in ("b_stress_vs_a_control", "a_stress_vs_b_control"):
            if "lfcSE" not in t[key].columns:
                raise KeyError("equivalence_margin requires an 'lfcSE' column in the matching contrasts")
        eq_b = equivalence_test(t["b_stress_vs_a_control"]["log2FoldChange"], t["b_stress_vs_a_control"]["lfcSE"],
                                opts.equivalence_margin)["padj"].to_numpy()
        eq_a = equivalence_test(t["a_stress_vs_b_control"]["log2FoldChange"], t["a_stress_vs_b_control"]["lfcSE"],
                                opts.equivalence_margin)["padj"].to_numpy()
        with np.errstate(invalid="ignore"):
            same_b = (~de_mb) & (eq_b < opts.alpha)
            same_a = (~de_ma) & (eq_a < opts.alpha)
    else:
        same_b = ~de_mb
        same_a = ~de_ma

    model = np.full(n, "", dtype=object)
    case = np.zeros(n, dtype=int)
    flags = [[] for _ in range(n)]

    only_b = ~de_a & de_b
    only_a = de_a & ~de_b
    both = de_a & de_b
    none = ~de_a & ~de_b

    # --- only group B responds: stress-ready 1/2 or unique 3/4 -----------------
    b_up = s_b == 1
    model[only_b & same_b] = "stress_ready"
    case[only_b & same_b] = np.where(b_up, 1, 2)[only_b & same_b]
    sel = only_b & de_mb
    model[sel] = "unique"
    case[sel] = np.where(b_up, 3, 4)[sel]
    sel = only_b & ~de_mb & ~same_b
    model[sel] = "unresolved"
    case[sel] = np.where(b_up, 1, 2)[sel]

    # --- only group A responds: stress-ready 3/4 or unique 1/2 -----------------
    a_up = s_a == 1
    model[only_a & same_a] = "stress_ready"
    case[only_a & same_a] = np.where(a_up, 3, 4)[only_a & same_a]
    sel = only_a & de_ma
    model[sel] = "unique"
    case[sel] = np.where(a_up, 1, 2)[sel]
    sel = only_a & ~de_ma & ~same_a
    model[sel] = "unresolved"
    case[sel] = np.where(a_up, 3, 4)[sel]

    # --- both respond: shared or opposite ----------------------------------------
    sel = both & (s_a == 1) & (s_b == 1)
    model[sel], case[sel] = "shared", 1
    sel = both & (s_a == -1) & (s_b == -1)
    model[sel], case[sel] = "shared", 2
    sel = both & (s_a == 1) & (s_b == -1)
    model[sel], case[sel] = "opposite", 1
    sel = both & (model == "")  # remaining sign combinations, as in the published script
    model[sel], case[sel] = "opposite", 2

    # --- no response ----------------------------------------------------------------
    sel = none & de_bas & (s_bas == 1)
    model[sel], case[sel] = "no_response", 1
    sel = none & de_bas & (s_bas == -1)
    model[sel], case[sel] = "no_response", 2
    sel = none & (model == "")
    model[sel], case[sel] = "no_response", 3

    # --- optional: require a significant basal difference for stress-ready -------
    if opts.require_basal_difference:
        sr = model == "stress_ready"
        expected = np.select([case == 1, case == 2, case == 3, case == 4], [-1, 1, 1, -1], 0)
        ok = de_bas & (s_bas == expected)
        demote = sr & ~ok
        model[demote] = "unique"
        case[demote] = np.select([case == 1, case == 2, case == 3, case == 4], [3, 4, 1, 2], 0)[demote]
        for i in np.flatnonzero(demote):
            flags[i].append("basal_difference_not_significant")

    # --- missing values ---------------------------------------------------------------
    consulted_na = na_a | na_b
    consulted_na |= only_b & na_mb
    consulted_na |= only_a & na_ma
    consulted_na |= none & na_bas
    if opts.require_basal_difference:
        consulted_na |= (only_a | only_b) & na_bas
    if opts.na_action == "unclassified":
        model[consulted_na] = "unclassified"
        case[consulted_na] = 0
    else:
        for i in np.flatnonzero(consulted_na):
            flags[i].append("na_treated_as_not_de")

    defs = model_definitions(group_names, include_extra=True).set_index(["model", "case"])
    keys = list(zip(model, case))
    out = pd.DataFrame(
        {
            "label": [m if m == "unclassified" else f"{m}_{c}" for m, c in keys],
            "model": model.astype(str),
            "case": case,
            "primed_group": [defs.loc[k, "primed_group"] for k in keys],
            "responding_group": [defs.loc[k, "responding_group"] for k in keys],
            "description": [defs.loc[k, "description"] for k in keys],
            "flags": [";".join(f) for f in flags],
        },
        index=index,
    )
    out.index.name = "feature"
    if keep_statistics:
        for key in CONTRAST_KEYS:
            out[f"{key}_log2FoldChange"] = pd.to_numeric(t[key]["log2FoldChange"], errors="coerce").to_numpy()
            out[f"{key}_padj"] = pd.to_numeric(t[key]["padj"], errors="coerce").to_numpy()
    out.attrs["options"] = opts
    out.attrs["group_names"] = tuple(group_names)
    return out


# ----------------------------------------------------------------------------
# Summaries
# ----------------------------------------------------------------------------

def summarize_classification(classification: pd.DataFrame) -> pd.DataFrame:
    """Counts per model and case.

    ``pct_of_responsive`` uses as denominator the features that respond to the
    stress in at least one group (stress-ready, shared, unique, opposite and
    unresolved); ``pct_of_all`` uses every classified feature. Reporting both
    avoids confusing the two denominators.
    """
    counts = classification.groupby(["model", "case"]).size().rename("n").reset_index()
    order = {m: i for i, m in enumerate(("stress_ready", "shared", "unique", "opposite", "no_response", "unresolved", "unclassified"))}
    counts["_o"] = counts["model"].map(order)
    counts = counts.sort_values(["_o", "case"]).drop(columns="_o").reset_index(drop=True)
    counts.insert(0, "label", [m if m == "unclassified" else f"{m}_{c}" for m, c in zip(counts["model"], counts["case"])])
    n_resp = int(classification["model"].isin(RESPONSIVE_MODELS).sum())
    n_all = len(classification)
    counts["pct_of_responsive"] = np.where(
        counts["model"].isin(RESPONSIVE_MODELS), 100 * counts["n"] / max(n_resp, 1), np.nan
    )
    counts["pct_of_all"] = 100 * counts["n"] / max(n_all, 1)
    return counts


def primed_asymmetry(classification: pd.DataFrame) -> pd.DataFrame:
    """Compare how often single-group responses end in a stress-ready state.

    Among features that respond to stress in only one group, the fraction that
    reaches the other group's control level is computed separately for A-only
    and B-only responders, and compared with Fisher's exact test. A larger raw
    number of stress-ready cases in one direction can simply reflect more
    responsive genes in that group; this test conditions on that.
    """
    m, c = classification["model"], classification["case"]
    a_sr = int(((m == "stress_ready") & c.isin([3, 4])).sum())
    a_un = int(((m == "unique") & c.isin([1, 2])).sum())
    b_sr = int(((m == "stress_ready") & c.isin([1, 2])).sum())
    b_un = int(((m == "unique") & c.isin([3, 4])).sum())
    names = classification.attrs.get("group_names", ("A", "B"))
    table = pd.DataFrame(
        {
            "responding_group": [names[0], names[1]],
            "primed_group": [names[1], names[0]],
            "stress_ready": [a_sr, b_sr],
            "unique": [a_un, b_un],
        }
    )
    table["fraction_stress_ready"] = table["stress_ready"] / (table["stress_ready"] + table["unique"]).replace(0, np.nan)
    odds, p = fisher_exact([[a_sr, a_un], [b_sr, b_un]])
    table.attrs["odds_ratio"] = odds
    table.attrs["fisher_pvalue"] = p
    table["odds_ratio"] = odds
    table["fisher_pvalue"] = p
    return table
