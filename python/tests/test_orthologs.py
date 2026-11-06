import numpy as np
import pandas as pd
import pytest

from altiro import best_hits, merge_ortholog_counts, reciprocal_best_hits


def _blast(rows):
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    return pd.DataFrame([dict(zip(cols, r)) for r in rows])


def test_rbh_basic(tmp_path):
    ab = _blast([
        ("a1.1", "b1.1", 90, 100, 1, 0, 1, 100, 1, 100, 1e-50, 300),
        ("a1.1", "b2.1", 60, 100, 1, 0, 1, 100, 1, 100, 1e-20, 120),
        ("a2.1", "b2.1", 80, 100, 1, 0, 1, 100, 1, 100, 1e-40, 250),
        ("a3.1", "b2.1", 70, 100, 1, 0, 1, 100, 1, 100, 1e-30, 200),
    ])
    ba = _blast([
        ("b1.1", "a1.1", 90, 100, 1, 0, 1, 100, 1, 100, 1e-50, 300),
        ("b2.1", "a2.1", 80, 100, 1, 0, 1, 100, 1, 100, 1e-40, 250),
        ("b2.1", "a3.1", 70, 100, 1, 0, 1, 100, 1, 100, 1e-30, 200),
    ])
    p_ab, p_ba = tmp_path / "ab.tsv", tmp_path / "ba.tsv"
    ab.to_csv(p_ab, sep="\t", header=False, index=False)
    ba.to_csv(p_ba, sep="\t", header=False, index=False)
    rbh = reciprocal_best_hits(p_ab, p_ba, strip_pattern=r"\.\d+$")
    assert list(zip(rbh["gene_a"], rbh["gene_b"])) == [("a1", "b1"), ("a2", "b2")]
    assert list(rbh["ortholog_id"]) == ["ortholog_1", "ortholog_2"]


def test_best_hits_tie_break():
    df = _blast([("q", "s2", 1, 1, 1, 1, 1, 1, 1, 1, 1e-5, 50), ("q", "s1", 1, 1, 1, 1, 1, 1, 1, 1, 1e-5, 50)])
    assert best_hits(df)["sseqid"].iloc[0] == "s1"


def test_merge_and_length_normalization():
    ca = pd.DataFrame({"a_c1": [100, 10], "a_s1": [200, 20]}, index=["a1", "a2"])
    cb = pd.DataFrame({"b_c1": [50, 5], "b_s1": [60, 6]}, index=["b1", "b2"])
    orth = pd.DataFrame({"ortholog_id": ["o1", "o2"], "gene_a": ["a1", "a2"], "gene_b": ["b1", "b2"]})
    m = merge_ortholog_counts(ca, cb, orth)
    assert m.shape == (2, 4) and m.loc["o1", "b_s1"] == 60
    m2 = merge_ortholog_counts(ca, cb, orth, lengths_a=pd.Series({"a1": 2000, "a2": 500}),
                               lengths_b=pd.Series({"b1": 1000, "b2": 4000}), length_normalize=True)
    assert m2.loc["o1", "a_c1"] == 50 and m2.loc["o2", "a_c1"] == 20 and m2.loc["o2", "b_c1"] == 1


def test_merge_rejects_non_one_to_one():
    ca = pd.DataFrame({"x": [1, 2]}, index=["a1", "a2"])
    cb = pd.DataFrame({"y": [1]}, index=["b1"])
    orth = pd.DataFrame({"ortholog_id": ["o1", "o2"], "gene_a": ["a1", "a2"], "gene_b": ["b1", "b1"]})
    with pytest.raises(ValueError):
        merge_ortholog_counts(ca, cb, orth)
