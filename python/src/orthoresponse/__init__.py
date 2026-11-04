"""orthoresponse: classify expression response models of orthologs (or genes) in 2 x 2 designs.

Two groups (species, genotypes, cell lines) x two conditions (control, stress).
Each feature is assigned to one of five models -- stress-ready, shared, unique,
opposite or no response -- with the decision tree of Contreras-Riquelme et al.
(2025) PLoS One 20(5): e0324724.
"""

__version__ = "0.2.0"

from .classify import (  # noqa: E402
    CONTRAST_KEYS,
    CONTRAST_LABELS,
    ClassificationOptions,
    adjust_bh,
    classification_options,
    classify_responses,
    equivalence_test,
    model_definitions,
    primed_asymmetry,
    summarize_classification,
)
from .datasets import ExampleData, example_files, export_example, list_examples, load_example  # noqa: E402
from .contrasts import ContrastSet, build_design, read_contrasts, run_contrasts, write_contrasts  # noqa: E402
from .orthologs import best_hits, merge_ortholog_counts, read_blast_tabular, reciprocal_best_hits  # noqa: E402
from .pipeline import OrthoResponseResult, run_orthoresponse  # noqa: E402

__all__ = [
    "__version__",
    "CONTRAST_KEYS",
    "CONTRAST_LABELS",
    "ClassificationOptions",
    "ContrastSet",
    "ExampleData",
    "OrthoResponseResult",
    "adjust_bh",
    "best_hits",
    "build_design",
    "classification_options",
    "classify_responses",
    "equivalence_test",
    "example_files",
    "export_example",
    "list_examples",
    "load_example",
    "merge_ortholog_counts",
    "model_definitions",
    "primed_asymmetry",
    "read_blast_tabular",
    "read_contrasts",
    "reciprocal_best_hits",
    "run_contrasts",
    "run_orthoresponse",
    "summarize_classification",
    "write_contrasts",
]
