from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from orthoresponse import (
    CONTRAST_KEYS,
    adjust_bh,
    classification_options,
    classify_responses,
    equivalence_test,
    model_definitions,
    primed_asymmetry,
    read_contrasts,
    summarize_classification,
)

from conftest import random_contrasts
from reference_notebook import PAPER_FILES, classify_like_notebook

DATA = Path(__file__).parent / "data"


def _write_paper_csvs(tables, folder):
    """Write tables in the exact layout produced by the published DESeq2 script."""
    for key, fname in PAPER_FILES.items():
        df = tables[key].reset_index().rename(columns={"feature": "Row.names"})
        df.to_csv(folder / fname)  # leading unnamed index column, like write.csv in R


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_paper_preset_matches_published_notebook(tmp_path, seed):
    tables = random_contrasts(n=3000, seed=seed)
    _write_paper_csvs(tables, tmp_path)
    expected = classify_like_notebook(tmp_path)
    paths = {k: tmp_path / f for k, f in PAPER_FILES.items()}
    got = classify_responses(read_contrasts(paths, id_col="Row.names"))
    assert set(expected.index) == set(got.index)
    assert (got.loc[expected.index, "label"] == expected).all()


def test_every_model_is_reachable(contrasts):
    got = classify_responses(contrasts)
    labels = set(got["label"])
    assert labels == set(model_definitions()["label"])


def test_single_gene_rules():
    def one(a, b, basal, mb, ma):
        tabs = {}
        for key, (padj, lfc) in zip(CONTRAST_KEYS, (a, b, basal, mb, ma)):
            tabs[key] = pd.DataFrame({"padj": [padj], "log2FoldChange": [lfc]}, index=["g"])
        return classify_responses(tabs)["label"].iloc[0]

    ns, de = 0.5, 0.001
    assert one((ns, 0.1), (de, 2), (de, -2), (ns, 0.1), (de, 2)) == "stress_ready_1"
    assert one((ns, 0.1), (de, -2), (de, 2), (ns, 0.1), (de, -2)) == "stress_ready_2"
    assert one((de, 2), (ns, 0.1), (de, 2), (de, -2), (ns, 0.1)) == "stress_ready_3"
    assert one((de, -2), (ns, 0.1), (de, -2), (de, 2), (ns, 0.1)) == "stress_ready_4"
    assert one((de, 2), (ns, 0.1), (ns, 0.1), (de, 2), (de, 2)) == "unique_1"
    assert one((ns, 0.1), (de, -2), (ns, 0.1), (de, -2), (ns, 0.1)) == "unique_4"
    assert one((de, 2), (de, 2), (ns, 0), (de, 2), (de, 2)) == "shared_1"
    assert one((de, -2), (de, 2), (ns, 0), (ns, 0), (ns, 0)) == "opposite_2"
    assert one((ns, 0), (ns, 0), (de, 1), (de, 1), (de, -1)) == "no_response_1"
    assert one((ns, 0), (ns, 0), (np.nan, 1), (ns, 0), (ns, 0)) == "no_response_3"


def test_strict_requires_basal_difference():
    tabs = {}
    # A induced, B flat, A stress ~ B control, but no basal difference -> unique_1 under strict
    vals = {"a_response": (0.001, 1.0), "b_response": (0.6, 0.1), "basal": (0.3, 0.4),
            "b_stress_vs_a_control": (0.2, 0.5), "a_stress_vs_b_control": (0.4, 0.3)}
    for k, (p, l) in vals.items():
        tabs[k] = pd.DataFrame({"padj": [p], "log2FoldChange": [l]}, index=["g"])
    assert classify_responses(tabs)["label"].iloc[0] == "stress_ready_3"
    strict = classify_responses(tabs, preset="strict")
    assert strict["label"].iloc[0] == "unique_1"
    assert "basal_difference_not_significant" in strict["flags"].iloc[0]


def test_na_handling(contrasts):
    paper = classify_responses(contrasts)
    strict = classify_responses(contrasts, na_action="unclassified")
    flagged = paper["flags"].str.contains("na_treated_as_not_de")
    assert flagged.sum() > 0
    assert (strict.loc[flagged, "model"] == "unclassified").all()
    assert (strict.loc[~flagged, "label"] == paper.loc[~flagged, "label"]).all()


def test_equivalence_only_moves_stress_ready(contrasts):
    paper = classify_responses(contrasts)
    eq = classify_responses(contrasts, equivalence_margin=0.5)
    changed = paper["label"] != eq["label"]
    assert (paper.loc[changed, "model"] == "stress_ready").all()
    assert (eq.loc[changed, "model"] == "unresolved").all()


def test_equivalence_test_against_closed_form():
    from scipy.stats import norm

    lfc, se, m = np.array([0.0, 0.3, 1.2]), np.array([0.1, 0.2, 0.2]), 0.5
    res = equivalence_test(lfc, se, m)
    expected = np.maximum(norm.sf((lfc + m) / se), norm.cdf((lfc - m) / se))
    assert np.allclose(res["pvalue"], expected)


def test_adjust_bh_matches_reference():
    p = np.array([0.01, 0.04, np.nan, 0.03, 0.2, 0.5])
    # values from R: p.adjust(c(0.01, 0.04, NA, 0.03, 0.2, 0.5), "BH")
    expected = np.array([0.05, 0.0666666667, np.nan, 0.0666666667, 0.25, 0.5])
    assert np.allclose(adjust_bh(p), expected, equal_nan=True)


def test_summary_denominators(contrasts):
    cls = classify_responses(contrasts)
    s = summarize_classification(cls)
    assert s["n"].sum() == len(cls)
    assert np.isclose(s.loc[s["pct_of_responsive"].notna(), "pct_of_responsive"].sum(), 100)
    assert np.isclose(s["pct_of_all"].sum(), 100)


def test_primed_asymmetry(contrasts):
    t = primed_asymmetry(classify_responses(contrasts))
    assert list(t.columns[:4]) == ["responding_group", "primed_group", "stress_ready", "unique"]
    assert 0 <= t["fisher_pvalue"].iloc[0] <= 1


def test_wide_table_input(contrasts):
    wide = pd.concat({k: v for k, v in contrasts.items()}, axis=1)
    wide.columns = [f"{k}_{c}" for k, c in wide.columns]
    a = classify_responses(contrasts)["label"]
    b = classify_responses(wide)["label"]
    assert (a == b).all()


def test_options_validation():
    with pytest.raises(ValueError):
        classification_options("nope")
    with pytest.raises(TypeError):
        classification_options("paper", not_an_option=1)
    with pytest.raises(ValueError):
        classification_options("paper", na_action="drop")


def test_shared_fixture_matches_expected():
    """Fixture shared with the R package: both implementations must give these labels."""
    paths = {k: DATA / f"contrast_{k}.tsv" for k in CONTRAST_KEYS}
    tables = read_contrasts(paths)
    expected = pd.read_csv(DATA / "expected_labels.tsv", sep="\t", index_col=0)
    for preset in ("paper", "strict"):
        got = classify_responses(tables, preset=preset)
        assert (got.loc[expected.index, "label"] == expected[preset]).all()
    got = classify_responses(tables, equivalence_margin=0.5)
    assert (got.loc[expected.index, "label"] == expected["paper_equivalence_0.5"]).all()
