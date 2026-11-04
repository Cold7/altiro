"""Command-line interface: ``orthoresponse <command> --help``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from . import __version__
from .classify import CONTRAST_KEYS, classify_responses, model_definitions, primed_asymmetry, summarize_classification


def _read_table(path: str, index_col=0) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(p, index_col=index_col)
    sep = "," if p.name.endswith((".csv", ".csv.gz")) else "\t"
    return pd.read_csv(p, sep=sep, index_col=index_col)


def _add_classification_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("classification")
    g.add_argument("--preset", default="paper", choices=["paper", "strict"])
    g.add_argument("--alpha", type=float, default=0.05)
    g.add_argument("--min-abs-lfc", type=float, default=None)
    g.add_argument("--na-action", choices=["not_de", "unclassified"], default=None)
    g.add_argument("--require-basal-difference", action="store_true", default=None)
    g.add_argument("--equivalence-margin", type=float, default=None)
    g.add_argument("--group-names", nargs=2, metavar=("A", "B"), default=None)


def _overrides(a) -> dict:
    return {
        "alpha": a.alpha,
        "min_abs_lfc": a.min_abs_lfc,
        "na_action": a.na_action,
        "require_basal_difference": a.require_basal_difference,
        "equivalence_margin": a.equivalence_margin,
    }


def _add_plot_args(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("plots")
    g.add_argument("--no-plots", action="store_true", help="do not save model_patterns / model_summary figures")
    g.add_argument("--condition-labels", nargs=2, metavar=("CONTROL", "STRESS"), default=("C", "S"),
                   help="tick labels for the two conditions in model_patterns (default: C S)")
    g.add_argument("--plot-formats", nargs="+", default=["png", "pdf"], help="figure formats (default: png pdf)")


def _write_outputs(cls: pd.DataFrame, outdir: Path, a=None) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    cls.to_csv(outdir / "classification.tsv", sep="\t")
    summ = summarize_classification(cls)
    summ.to_csv(outdir / "summary.tsv", sep="\t", index=False)
    primed_asymmetry(cls).to_csv(outdir / "primed_asymmetry.tsv", sep="\t", index=False)
    if a is not None and not a.no_plots:
        from .pipeline import save_plots

        save_plots(cls, outdir, condition_labels=a.condition_labels, formats=a.plot_formats)
    with pd.option_context("display.width", 120, "display.max_columns", 10):
        print(summ.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print(f"\nResults written to {outdir}/")


def cmd_run(a) -> None:
    from .pipeline import run_orthoresponse

    counts = _read_table(a.counts)
    samples = _read_table(a.samples)
    res = run_orthoresponse(
        counts, samples, group_col=a.group_col, condition_col=a.condition_col, group_a=a.group_a,
        group_b=a.group_b, control=a.control, stress=a.stress, preset=a.preset,
        min_total_count=a.min_total_count, filter_alpha=a.filter_alpha, round_counts=a.round_counts,
        **_overrides(a),
    )
    res.save(a.out, plots=not a.no_plots, condition_labels=a.condition_labels, formats=a.plot_formats)
    _write_outputs(res.classification, Path(a.out))


def cmd_classify(a) -> None:
    from .contrasts import read_contrasts

    paths = {k: getattr(a, k) for k in CONTRAST_KEYS}
    tables = read_contrasts(paths, id_col=a.id_col)
    names = tuple(a.group_names) if a.group_names else ("A", "B")
    cls = classify_responses(tables, preset=a.preset, group_names=names, **_overrides(a))
    _write_outputs(cls, Path(a.out), a)


def cmd_rbh(a) -> None:
    from .orthologs import reciprocal_best_hits

    rbh = reciprocal_best_hits(a.ab, a.ba, max_evalue=a.max_evalue, min_bitscore=a.min_bitscore,
                               strip_pattern=a.strip_pattern)
    rbh.to_csv(a.out, sep="\t", index=False)
    print(f"{len(rbh)} reciprocal best hits written to {a.out}")


def cmd_merge(a) -> None:
    from .orthologs import merge_ortholog_counts

    ca, cb = _read_table(a.counts_a), _read_table(a.counts_b)
    orth = pd.read_csv(a.orthologs, sep="\t")
    la = lb = None
    if a.lengths_a and a.lengths_b:
        la = _read_table(a.lengths_a).iloc[:, 0]
        lb = _read_table(a.lengths_b).iloc[:, 0]
    merged = merge_ortholog_counts(ca, cb, orth, lengths_a=la, lengths_b=lb, length_normalize=a.length_normalize)
    merged.to_csv(a.out, sep="\t")
    print(f"Merged matrix with {merged.shape[0]} orthologs x {merged.shape[1]} samples written to {a.out}")


def cmd_example(a) -> None:
    from .datasets import _EXAMPLES, export_example

    out = export_example(a.out, a.name)
    p = _EXAMPLES[a.name]["params"]
    print(f"Wrote {out['counts']} and {out['samples']}")
    print("Run it with:")
    print(f"  orthoresponse run --counts {out['counts']} --samples {out['samples']} "
          f"--group-col {p['group_col']} --condition-col {p['condition_col']} --group-a {p['group_a']} "
          f"--group-b {p['group_b']} --control {p['control']} --stress {p['stress']} --out results")


def cmd_models(a) -> None:
    names = tuple(a.group_names) if a.group_names else ("A", "B")
    print(model_definitions(names, include_extra=a.all).to_string(index=False))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="orthoresponse", description="Classify expression response models in 2 x 2 designs.")
    p.add_argument("--version", action="version", version=f"orthoresponse {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    r = sub.add_parser("run", help="counts + sample table -> PyDESeq2 -> classification")
    r.add_argument("--counts", required=True, help="raw counts (features x samples; TSV/CSV/XLSX, first column = id)")
    r.add_argument("--samples", required=True, help="sample table (first column = sample name)")
    r.add_argument("--group-col", required=True)
    r.add_argument("--condition-col", required=True)
    r.add_argument("--group-a", required=True, help="reference group (e.g. sensitive)")
    r.add_argument("--group-b", required=True, help="comparison group (e.g. tolerant)")
    r.add_argument("--control", required=True)
    r.add_argument("--stress", required=True)
    r.add_argument("--min-total-count", type=float, default=10)
    r.add_argument("--filter-alpha", type=float, default=None)
    r.add_argument("--round-counts", action="store_true")
    r.add_argument("--out", required=True, help="output directory")
    _add_classification_args(r)
    _add_plot_args(r)
    r.set_defaults(func=cmd_run)

    c = sub.add_parser("classify", help="classify from five precomputed contrast tables")
    for k in CONTRAST_KEYS:
        c.add_argument(f"--{k.replace('_', '-')}", dest=k, required=True, help=f"table for {k}")
    c.add_argument("--id-col", default=None, help="feature id column (default: first column)")
    c.add_argument("--out", required=True)
    _add_classification_args(c)
    _add_plot_args(c)
    c.set_defaults(func=cmd_classify)

    b = sub.add_parser("rbh", help="reciprocal best hits from two BLAST -outfmt 6 files")
    b.add_argument("--ab", required=True)
    b.add_argument("--ba", required=True)
    b.add_argument("--max-evalue", type=float, default=None)
    b.add_argument("--min-bitscore", type=float, default=None)
    b.add_argument("--strip-pattern", default=None)
    b.add_argument("--out", required=True)
    b.set_defaults(func=cmd_rbh)

    m = sub.add_parser("merge", help="merge two species' counts through an ortholog table")
    m.add_argument("--counts-a", required=True)
    m.add_argument("--counts-b", required=True)
    m.add_argument("--orthologs", required=True)
    m.add_argument("--lengths-a", default=None)
    m.add_argument("--lengths-b", default=None)
    m.add_argument("--length-normalize", action="store_true")
    m.add_argument("--out", required=True)
    m.set_defaults(func=cmd_merge)

    e = sub.add_parser("example", help="write the bundled example data (counts + samples) to a folder")
    e.add_argument("--name", default="rice_salt", choices=["rice_salt", "ovarian_brequinar"])
    e.add_argument("--out", default="orthoresponse_example")
    e.set_defaults(func=cmd_example)

    d = sub.add_parser("models", help="print the definition of every model and case")
    d.add_argument("--group-names", nargs=2, metavar=("A", "B"), default=None)
    d.add_argument("--all", action="store_true", help="include unresolved/unclassified")
    d.set_defaults(func=cmd_models)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
