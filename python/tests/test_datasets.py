import pandas as pd

import orthoresponse as orr


def test_load_rice_example():
    ex = orr.load_example("rice_salt")
    assert ex.counts.shape == (36844, 12)
    assert list(ex.counts.columns) == list(ex.samples.index)
    assert int(ex.counts.to_numpy().sum()) == 198199444
    d = orr.build_design(ex.samples, ex.params["group_col"], ex.params["condition_col"], ex.params["group_a"],
                         ex.params["group_b"], ex.params["control"], ex.params["stress"])
    assert d["cell"].value_counts().to_dict() == {"A_control": 3, "A_stress": 3, "B_control": 3, "B_stress": 3}
    assert "rice_salt" in orr.list_examples()["name"].tolist()


def test_export_example(tmp_path):
    out = orr.export_example(tmp_path / "ex")
    counts = pd.read_csv(out["counts"], sep="\t", index_col=0)
    samples = pd.read_csv(out["samples"], sep="\t", index_col=0)
    assert counts.shape == (36844, 12) and len(samples) == 12


def test_cli_example(tmp_path, capsys):
    from orthoresponse.cli import main
    assert main(["example", "--out", str(tmp_path / "ex")]) == 0
    assert "orthoresponse run --counts" in capsys.readouterr().out


def test_load_ovarian_example():
    ex = orr.load_example("ovarian_brequinar")
    assert ex.counts.shape == (59453, 12)
    assert list(ex.counts.columns) == list(ex.samples.index)
    assert int(ex.counts.to_numpy().sum()) == 281005562
    d = orr.build_design(ex.samples, ex.params["group_col"], ex.params["condition_col"], ex.params["group_a"],
                         ex.params["group_b"], ex.params["control"], ex.params["stress"])
    assert d["cell"].value_counts().to_dict() == {"A_control": 3, "A_stress": 3, "B_control": 3, "B_stress": 3}
    assert set(orr.list_examples()["name"]) == {"rice_salt", "ovarian_brequinar"}
