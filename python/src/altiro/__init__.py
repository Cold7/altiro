"""ALTIRO (Analysis of Latent Transcriptional Induction and Response Outcomes).

Classify expression response models of genes or orthologs in 2 x 2 designs.

Two groups (species, genotypes, cell lines) x two conditions (control, stress).
Each feature is assigned to one of five models -- primed (called stress-ready in
the original publication), shared, unique,
opposite or no response -- with the decision tree of Contreras-Riquelme et al.
(2025) PLoS One 20(5): e0324724.
"""

__version__ = "0.6.0"

from .classify import (  # noqa: E402
    CONTRAST_KEYS,
    CONTRAST_LABELS,
    PRIMED_MODELS,
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
from .pipeline import AltiroResult, format_asymmetry, run_altiro, save_plots  # noqa: E402

__all__ = [
    "__version__",
    "CONTRAST_KEYS",
    "CONTRAST_LABELS",
    "ClassificationOptions",
    "PRIMED_MODELS",
    "ContrastSet",
    "ExampleData",
    "AltiroResult",
    "adjust_bh",
    "best_hits",
    "build_design",
    "classification_options",
    "classify_responses",
    "equivalence_test",
    "example_files",
    "format_asymmetry",
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
    "run_altiro",
    "save_plots",
    "summarize_classification",
    "write_contrasts",
]
