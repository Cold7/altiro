import numpy as np
import pandas as pd
import pytest

from orthoresponse import CONTRAST_KEYS


def random_contrasts(n=2000, seed=0, na_frac=0.05):
    """Random contrast tables that visit every branch of the decision tree."""
    rng = np.random.default_rng(seed)
    idx = [f"ortholog_{i}" for i in range(n)]
    tables = {}
    for key in CONTRAST_KEYS:
        padj = rng.choice([0.001, 0.01, 0.049, 0.05, 0.2, 0.9], size=n)
        lfc = rng.normal(0, 1.5, size=n)
        lfc[rng.random(n) < 0.03] = 0.0           # exact zero fold changes
        lfc[rng.random(n) < 0.03] = rng.choice([0.004, -0.004, 0.006, -0.006], size=1)[0]
        padj[rng.random(n) < na_frac] = np.nan     # independent filtering / outliers
        se = rng.uniform(0.05, 1.0, size=n)
        tables[key] = pd.DataFrame({"baseMean": rng.uniform(5, 5000, n), "log2FoldChange": lfc,
                                    "lfcSE": se, "padj": padj}, index=pd.Index(idx, name="feature"))
    return tables


@pytest.fixture
def contrasts():
    return random_contrasts()
