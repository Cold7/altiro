#' Run the full workflow: DESeq2 contrasts and classification
#'
#' Convenience wrapper around [run_contrasts()], [classify_responses()],
#' [summarize_classification()] and [primed_asymmetry()].
#'
#' @inheritParams run_contrasts
#' @param preset `"paper"` (default) or `"strict"`.
#' @param ... Classification options passed to [classify_responses()], e.g.
#'   `equivalence_margin = 0.5`.
#' @return An object of class `orthoresponse` with elements `contrasts`,
#'   `classification`, `summary` and `asymmetry`.
#' @examples
#' \dontrun{
#' res <- orthoresponse(counts, samples, group_col = "line", condition_col = "treatment",
#'                      group_a = "PL12", group_b = "IL22", control = "Control", stress = "Salt")
#' res
#' save_orthoresponse(res, "results")
#' }
#' @export
orthoresponse <- function(counts, samples, group_col, condition_col, group_a, group_b, control, stress,
                          preset = c("paper", "strict"), min_total_count = 10, alpha = 0.05,
                          filter_alpha = NULL, round_counts = FALSE, quiet = TRUE, ...) {
  preset <- match.arg(preset)
  cs <- run_contrasts(counts, samples, group_col = group_col, condition_col = condition_col,
                      group_a = group_a, group_b = group_b, control = control, stress = stress,
                      min_total_count = min_total_count, alpha = alpha, filter_alpha = filter_alpha,
                      round_counts = round_counts, quiet = quiet)
  cls <- classify_responses(cs, preset = preset, group_names = cs$group_names, alpha = alpha, ...)
  structure(list(contrasts = cs, classification = cls, summary = summarize_classification(cls),
                 asymmetry = primed_asymmetry(cls)), class = "orthoresponse")
}

#' @export
print.orthoresponse <- function(x, ...) {
  cat("orthoresponse result: A =", x$contrasts$group_names[1], "| B =", x$contrasts$group_names[2], "\n")
  cat(nrow(x$classification), "features classified\n\n")
  s <- x$summary
  s$pct_of_responsive <- round(s$pct_of_responsive, 1)
  s$pct_of_all <- round(s$pct_of_all, 1)
  print(s, row.names = FALSE)
  a <- x$asymmetry
  cat(sprintf("\nStress-ready among single-group responders: %s-only %.1f%%, %s-only %.1f%% (Fisher p = %.3g)\n",
              a$responding_group[1], 100 * a$fraction_stress_ready[1],
              a$responding_group[2], 100 * a$fraction_stress_ready[2], a$fisher_pvalue[1]))
  invisible(x)
}

#' Save an orthoresponse result to disk
#'
#' Writes `classification.tsv`, `summary.tsv`, `primed_asymmetry.tsv`,
#' `normalized_counts.tsv` and the five contrast tables.
#'
#' @param x Output of [orthoresponse()].
#' @param outdir Output directory.
#' @return Invisibly, `outdir`.
#' @export
save_orthoresponse <- function(x, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  w <- function(d, f) utils::write.table(d, file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
  w(x$classification, "classification.tsv")
  w(x$summary, "summary.tsv")
  w(x$asymmetry, "primed_asymmetry.tsv")
  nc <- x$contrasts$normalized_counts
  w(data.frame(feature = rownames(nc), nc, check.names = FALSE), "normalized_counts.tsv")
  write_contrasts(x$contrasts, file.path(outdir, "contrasts"))
  invisible(outdir)
}
