#' Run the full workflow: DESeq2 contrasts and classification
#'
#' Convenience wrapper around [run_contrasts()], [classify_responses()],
#' [summarize_classification()] and [primed_asymmetry()].
#'
#' @inheritParams run_contrasts
#' @param preset `"paper"` (default) or `"strict"`.
#' @param group_names Display names for groups A and B in descriptions,
#'   summaries and figures (default: `group_a` and `group_b`, i.e. the values of
#'   the sample table).
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
                          filter_alpha = NULL, round_counts = FALSE, quiet = TRUE, group_names = NULL, ...) {
  preset <- match.arg(preset)
  cs <- run_contrasts(counts, samples, group_col = group_col, condition_col = condition_col,
                      group_a = group_a, group_b = group_b, control = control, stress = stress,
                      min_total_count = min_total_count, alpha = alpha, filter_alpha = filter_alpha,
                      round_counts = round_counts, quiet = quiet)
  if (!is.null(group_names)) {
    if (length(group_names) != 2) stop("group_names must have exactly two elements: c(A, B)")
    cs$group_names <- as.character(group_names)
  }
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
  cat(sprintf("\nSingle-group responses ending primed: %s-only %.1f%%, %s-only %.1f%% (Fisher p = %.3g)\n",
              a$responding_group[1], 100 * a$fraction_primed[1],
              a$responding_group[2], 100 * a$fraction_primed[2], a$fisher_pvalue[1]))
  invisible(x)
}

#' Save an orthoresponse result to disk
#'
#' Writes `classification.tsv`, `summary.tsv`, `primed_asymmetry.tsv`,
#' `normalized_counts.tsv`, the five contrast tables and, with `plots = TRUE`
#' (default), the figures `model_patterns.<fmt>` (the 15 patterns with counts)
#' and `model_summary.<fmt>` (bar chart).
#'
#' @param x Output of [orthoresponse()].
#' @param outdir Output directory.
#' @param plots Save the figures (default `TRUE`).
#' @param ... Passed to [save_plots()]: `condition_labels`, `model_labels`,
#'   `group_names`, `title`, `formats`, `width`, `height`, `res`.
#' @return Invisibly, `outdir`.
#' @examples
#' \dontrun{
#' save_orthoresponse(res, "results", condition_labels = c("Ctrl", "NaCl"))
#' }
#' @export
save_orthoresponse <- function(x, outdir, plots = TRUE, ...) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  w <- function(d, f) utils::write.table(d, file.path(outdir, f), sep = "\t", quote = FALSE, row.names = FALSE)
  w(x$classification, "classification.tsv")
  w(x$summary, "summary.tsv")
  w(x$asymmetry, "primed_asymmetry.tsv")
  nc <- x$contrasts$normalized_counts
  w(data.frame(feature = rownames(nc), nc, check.names = FALSE), "normalized_counts.tsv")
  write_contrasts(x$contrasts, file.path(outdir, "contrasts"))
  if (plots) save_plots(x$classification, outdir, ...)
  invisible(outdir)
}

#' Save the model figures of a classification
#'
#' Writes `model_patterns.<fmt>` ([plot_model_patterns()]) and
#' `model_summary.<fmt>` ([plot_summary()]) for each format in `formats`.
#'
#' @param classification Output of [classify_responses()].
#' @param outdir Output directory.
#' @param condition_labels,model_labels,group_names,title Passed to
#'   [plot_model_patterns()].
#' @param formats Any of `"png"` and `"pdf"` (default both).
#' @param width,height Size of `model_patterns` in inches.
#' @param res Resolution (dpi) of the PNG files.
#' @return Invisibly, the paths of the written files.
#' @export
save_plots <- function(classification, outdir, condition_labels = c("C", "S"), model_labels = NULL,
                       group_names = NULL, title = NULL, formats = c("png", "pdf"),
                       width = 10, height = 9.5, res = 150) {
  formats <- match.arg(formats, c("png", "pdf"), several.ok = TRUE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  n_bars <- nrow(summarize_classification(classification))
  figs <- list(
    model_patterns = list(w = width, h = height, f = function() {
      plot_model_patterns(classification, group_names = group_names, condition_labels = condition_labels,
                          model_labels = model_labels, title = title)
    }),
    model_summary = list(w = 6, h = 0.3 * n_bars + 1.2, f = function() plot_summary(classification))
  )
  paths <- character()
  for (name in names(figs)) {
    fg <- figs[[name]]
    for (fmt in formats) {
      path <- file.path(outdir, paste0(name, ".", fmt))
      if (fmt == "png") {
        grDevices::png(path, width = fg$w, height = fg$h, units = "in", res = res)
      } else {
        grDevices::pdf(path, width = fg$w, height = fg$h)
      }
      tryCatch(fg$f(), finally = grDevices::dev.off())
      paths <- c(paths, path)
    }
  }
  invisible(paths)
}
