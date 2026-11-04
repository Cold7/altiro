"""Example data shipped with the package.

``rice_salt``: rice leaves, GEO GSE167342. Salt-sensitive *japonica* cultivar
PL12 (group A) and the salt-tolerant introgression line IL22 (group B, derived
from FL478 x PL12, carries the *Saltol* QTL), control vs 24 h 80 mM NaCl,
3 replicates per group. Raw gene-level counts (36,844 genes x 12 samples) as
deposited by the authors in GEO (``GSE167342_Raw_counts.xlsx``); columns
S40-S72 are the authors' library names (see the sample table).

``ovarian_brequinar``: human ovarian cancer cells, GEO GSE283381. OVCAR5
cisplatin-sensitive cells (CisS, group A) and their cisplatin-resistant
derivative (CisR, group B), treated with vehicle or the DHODH inhibitor
brequinar (de novo pyrimidine synthesis), 3 replicates per group. Raw HTSeq
counts deposited by the authors (one file per sample in ``GSE283381_RAW.tar``),
merged into one matrix (59,453 genes, Ensembl GRCh38 release 91 gene names x
12 samples); columns are the authors' library names.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Dict, Union

import pandas as pd

__all__ = ["ExampleData", "list_examples", "load_example", "example_files", "export_example"]

_EXAMPLES = {
    "rice_salt": {
        "counts": "rice_GSE167342_counts.tsv.gz",
        "samples": "rice_GSE167342_samples.tsv",
        "params": {
            "group_col": "line",
            "condition_col": "treatment",
            "group_a": "PL12",
            "group_b": "IL22",
            "control": "Control",
            "stress": "Salt",
        },
        "description": (
            "Rice leaves (GEO GSE167342): salt-sensitive PL12 (A) vs salt-tolerant Saltol introgression "
            "line IL22 (B), control vs 24 h 80 mM NaCl, 3 replicates. Raw counts, 36,844 genes x 12 samples."
        ),
    },
    "ovarian_brequinar": {
        "counts": "ovarian_GSE283381_counts.tsv.gz",
        "samples": "ovarian_GSE283381_samples.tsv",
        "params": {
            "group_col": "line",
            "condition_col": "treatment",
            "group_a": "CisS",
            "group_b": "CisR",
            "control": "Vehicle",
            "stress": "Brequinar",
        },
        "description": (
            "Human ovarian cancer cells (GEO GSE283381): OVCAR5 cisplatin-sensitive (A) vs cisplatin-resistant "
            "(B), vehicle vs the DHODH inhibitor brequinar, 3 replicates. Raw HTSeq counts, 59,453 genes x 12 samples."
        ),
    },
}


@dataclass
class ExampleData:
    """An example data set: raw counts, sample table and the matching arguments.

    ``params`` holds the arguments for :func:`~orthoresponse.run_orthoresponse`,
    so the whole example runs with ``run_orthoresponse(ex.counts, ex.samples, **ex.params)``.
    """

    name: str
    counts: pd.DataFrame
    samples: pd.DataFrame
    params: Dict[str, str] = field(default_factory=dict)
    description: str = ""

    def __repr__(self) -> str:
        return (f"ExampleData({self.name!r}: counts {self.counts.shape[0]} x {self.counts.shape[1]}, "
                f"{len(self.samples)} samples; params={self.params})")


def list_examples() -> pd.DataFrame:
    """Names and descriptions of the bundled example data sets."""
    return pd.DataFrame([{"name": k, "description": v["description"]} for k, v in _EXAMPLES.items()])


def _get(name: str) -> dict:
    if name not in _EXAMPLES:
        raise ValueError(f"Unknown example {name!r}; available: {list(_EXAMPLES)}")
    return _EXAMPLES[name]


def example_files(name: str = "rice_salt") -> Dict[str, Path]:
    """Paths of the bundled files (``counts`` and ``samples``) of an example."""
    info = _get(name)
    base = resources.files("orthoresponse") / "data"
    return {k: Path(str(base / info[k])) for k in ("counts", "samples")}


def load_example(name: str = "rice_salt") -> ExampleData:
    """Load a bundled example data set.

    Examples
    --------
    >>> import orthoresponse as orr
    >>> ex = orr.load_example("rice_salt")
    >>> res = orr.run_orthoresponse(ex.counts, ex.samples, **ex.params)
    >>> print(res.summary)
    """
    info = _get(name)
    files = example_files(name)
    counts = pd.read_csv(files["counts"], sep="\t", index_col=0)
    samples = pd.read_csv(files["samples"], sep="\t", index_col=0)
    return ExampleData(name=name, counts=counts, samples=samples, params=dict(info["params"]),
                       description=info["description"])


def export_example(outdir: Union[str, Path] = ".", name: str = "rice_salt") -> Dict[str, Path]:
    """Copy the files of an example to ``outdir`` (counts as ``counts.tsv.gz``, plus ``samples.tsv``).

    Useful as a template for preparing your own input files.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    src = example_files(name)
    out = {"counts": outdir / "counts.tsv.gz", "samples": outdir / "samples.tsv"}
    for key in out:
        shutil.copyfile(src[key], out[key])
    return out
