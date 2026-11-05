import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pydeseq2")

from orthoresponse import build_design, run_orthoresponse  # noqa: E402
from orthoresponse.cli import main as cli_main  # noqa: E402


def _simulate(seed=3):
    rng = np.random.default_rng(seed)
    n = 1500
    base = np.exp(rng.uniform(np.log(30), np.log(3000), n))
    # columns: A_control, A_stress, B_control, B_stress
    f = np.ones((n, 4))
    f[:150, [1, 2, 3]] = 4.0        # primed_3: A induced to B's constitutive level
    f[150:300, 1] = 4.0             # unique_1: induced only in A
    f[300:450, [1, 3]] = 0.25       # shared_2
    f[450:600, [2, 3]] = 3.0        # no_response_1
    f[600:675, 3] = 4.0             # unique_3
    f[675:750, [1, 3]] = 4.0        # shared_1 (balances fold-change directions)
    mu = np.repeat(base[:, None] * f, 3, axis=1)
    disp = 0.03
    counts = rng.negative_binomial(1 / disp, 1 / (1 + disp * mu))
    samples = [f"{g}_{c}_{r}" for g in ("sens", "tol") for c in ("ctrl", "salt") for r in (1, 2, 3)]
    counts = pd.DataFrame(counts, index=[f"gene{i}" for i in range(n)], columns=samples)
    meta = pd.DataFrame({"line": [s.split("_")[0] for s in samples],
                         "treatment": [s.split("_")[1] for s in samples]}, index=samples)
    return counts, meta


def test_design_mapping():
    _, meta = _simulate()
    d = build_design(meta, "line", "treatment", "sens", "tol", "ctrl", "salt")
    assert d["cell"].value_counts().to_dict() == {"A_control": 3, "A_stress": 3, "B_control": 3, "B_stress": 3}


def test_pipeline_recovers_simulated_models(tmp_path):
    counts, meta = _simulate()
    res = run_orthoresponse(counts, meta, group_col="line", condition_col="treatment", group_a="sens",
                            group_b="tol", control="ctrl", stress="salt")
    lab = res.classification["label"]
    renamed = run_orthoresponse(counts.iloc[:300], meta, group_col="line", condition_col="treatment", group_a="sens",
                                group_b="tol", control="ctrl", stress="salt", group_names=("Sensitive", "Tolerant"))
    assert str(renamed).startswith("orthoresponse result: A = Sensitive | B = Tolerant")
    assert set(renamed.classification["primed_group"]) <= {"Sensitive", "Tolerant", ""}
    assert list(renamed.asymmetry["responding_group"]) == ["Sensitive", "Tolerant"]
    assert (lab.iloc[:150] == "primed_3").mean() > 0.9
    assert (lab.iloc[150:300] == "unique_1").mean() > 0.9
    assert (lab.iloc[300:450] == "shared_2").mean() > 0.9
    assert (lab.iloc[450:600] == "no_response_1").mean() > 0.9
    text = str(res)
    assert text.startswith("orthoresponse result: A = sens | B = tol")
    assert "Single-group responses ending primed: sens-only" in text and "(Fisher p = " in text
    assert repr(res) == text
    out = res.save(tmp_path / "res", condition_labels=("Ctrl", "Salt"))
    assert (out / "classification.tsv").exists() and (out / "contrasts" / "contrast_basal.tsv").exists()
    for f in ("model_patterns.png", "model_patterns.pdf", "model_summary.png", "model_summary.pdf"):
        assert (out / f).stat().st_size > 1000, f
    out2 = res.save(tmp_path / "res2", plots=False)
    assert not (out2 / "model_patterns.png").exists()


def test_cli_run(tmp_path):
    counts, meta = _simulate()
    counts.to_csv(tmp_path / "counts.tsv", sep="\t")
    meta.to_csv(tmp_path / "samples.tsv", sep="\t")
    rc = cli_main(["run", "--counts", str(tmp_path / "counts.tsv"), "--samples", str(tmp_path / "samples.tsv"),
                   "--group-col", "line", "--condition-col", "treatment", "--group-a", "sens", "--group-b", "tol",
                   "--control", "ctrl", "--stress", "salt", "--group-names", "Sensitive", "Tolerant",
                   "--out", str(tmp_path / "out")])
    assert rc == 0
    cls = pd.read_csv(tmp_path / "out" / "classification.tsv", sep="\t", index_col=0)
    assert "label" in cls.columns and len(cls) > 1000
    assert "Tolerant" in set(cls["primed_group"])
    assert (tmp_path / "out" / "model_patterns.png").exists() and (tmp_path / "out" / "model_summary.pdf").exists()
