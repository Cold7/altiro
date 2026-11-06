#' altiro: Analysis of Latent Transcriptional Induction and Response Outcomes (ALTIRO)
#'
#' Classifies orthologs (cross-species) or genes (within species) into
#' primed (called stress-ready in the original publication), shared, unique,
#' opposite and no-response models from five
#' differential-expression contrasts between two groups (species, genotypes,
#' cell lines) under control and stress conditions. The default rules
#' reproduce the classification rules of Contreras-Riquelme et al. (2025)
#' \doi{10.1371/journal.pone.0324724}; optional stricter rules handle missing
#' statistics, require a basal difference, or test equivalence explicitly.
#'
#' Main functions:
#' * [altiro()] - counts to classification in one call.
#' * [run_contrasts()] - DESeq2 fit and the five contrasts.
#' * [classify_responses()] - the classification rules.
#' * [summarize_classification()], [primed_asymmetry()] - summaries.
#' * [reciprocal_best_hits()], [merge_ortholog_counts()] - cross-species input.
#' * [plot_model_patterns()], [plot_summary()] - figures.
#'
#' @keywords internal
"_PACKAGE"
