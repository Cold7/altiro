"""Regenerate the cross-language test fixture.

Writes five contrast tables and the expected labels (computed with the Python
implementation, whose ``paper`` preset is itself validated against the
published notebook) into ``python/tests/data`` and ``r-package/inst/extdata``.
The R test suite checks that it reproduces the same labels.

Usage: python tools/make_fixture.py
"""

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python" / "src"))
sys.path.insert(0, str(ROOT / "python" / "tests"))

import pandas as pd  # noqa: E402

from conftest import random_contrasts  # noqa: E402
from orthoresponse import CONTRAST_KEYS, classify_responses  # noqa: E402

PY = ROOT / "python" / "tests" / "data"
R = ROOT / "r-package" / "inst" / "extdata"


def main():
    PY.mkdir(parents=True, exist_ok=True)
    R.mkdir(parents=True, exist_ok=True)
    tables = random_contrasts(n=800, seed=42, na_frac=0.06)
    for key in CONTRAST_KEYS:
        tables[key].to_csv(PY / f"contrast_{key}.tsv", sep="\t", float_format="%.10g")
    # re-read so that the labels refer to exactly the written values
    from orthoresponse import read_contrasts

    t = read_contrasts({k: PY / f"contrast_{k}.tsv" for k in CONTRAST_KEYS})
    exp = pd.DataFrame({
        "paper": classify_responses(t)["label"],
        "strict": classify_responses(t, preset="strict")["label"],
        "paper_equivalence_0.5": classify_responses(t, equivalence_margin=0.5)["label"],
    })
    exp.to_csv(PY / "expected_labels.tsv", sep="\t")
    for f in PY.glob("*.tsv"):
        shutil.copy(f, R / f.name)
    print(exp.apply(lambda c: c.value_counts()).fillna(0).astype(int).to_string())


if __name__ == "__main__":
    main()
