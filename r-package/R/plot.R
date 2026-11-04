.patterns <- list(
  stress_ready = list(c(0, 0, -1, 0), c(0, 0, 1, 0), c(0, 1, 1, 1), c(0, -1, -1, -1)),
  shared = list(c(0, 1, 0, 1), c(0, -1, 0, -1)),
  unique = list(c(0, 1, 0, 0), c(0, -1, 0, 0), c(0, 0, 0, 1), c(0, 0, 0, -1)),
  opposite = list(c(0, 1, 0, -1), c(0, -1, 0, 1)),
  no_response = list(c(0, 0, 1, 1), c(0, 0, -1, -1), c(0, 0, 0, 0))
)
# panel column (1-4) of each case; shared and opposite use the outer columns, as in the original Fig. 3
.columns <- list(stress_ready = 1:4, shared = c(1, 4), unique = 1:4, opposite = c(1, 4), no_response = 1:3)
.default_model_labels <- c(stress_ready = "Stress-ready", shared = "Shared response", unique = "Unique response",
                           opposite = "Opposite response", no_response = "No response")
.col_a <- "#1f77b4"
.col_b <- "#ff7f0e"

#' Plot the expression model patterns
#'
#' Draws the 15 schematic patterns (left half of each panel group A, right
#' half group B; within each half, control on the left and stress on the
#' right), annotated with counts and percentages of responsive features when a
#' classification is supplied. Uses base graphics.
#'
#' @param classification Optional output of [classify_responses()].
#' @param group_names Names for the legend; taken from `classification` when
#'   available.
#' @param condition_labels Character vector of length 2, `c(control, stress)`,
#'   used as tick labels above every panel. Default `c("C", "S")`; e.g.
#'   `c("C", "D")` for drought as in the original publication, or
#'   `c("Ctrl", "NaCl")`.
#' @param model_labels Row labels: a named character vector (or list) with
#'   any of the names `stress_ready`, `shared`, `unique`, `opposite`,
#'   `no_response`. Only the names you give are changed, e.g.
#'   `c(stress_ready = "Primed", no_response = "Unchanged")`.
#' @param title Figure title. By default it explains the condition labels.
#' @return Invisibly, `NULL`.
#' @examples
#' plot_model_patterns()
#' plot_model_patterns(condition_labels = c("Ctrl", "Drought"), group_names = c("Sly", "Spenn"),
#'                     model_labels = c(stress_ready = "Primed", no_response = "Unchanged"))
#' @export
plot_model_patterns <- function(classification = NULL, group_names = NULL, condition_labels = c("C", "S"),
                                model_labels = NULL, title = NULL) {
  row_labels <- .default_model_labels
  if (!is.null(model_labels)) {
    model_labels <- unlist(model_labels)
    bad <- setdiff(names(model_labels), names(row_labels))
    if (is.null(names(model_labels)) || any(names(model_labels) == "") || length(bad)) {
      stop("model_labels must be named with any of: ", paste(names(row_labels), collapse = ", "))
    }
    row_labels[names(model_labels)] <- as.character(model_labels)
  }
  condition_labels <- as.character(condition_labels)
  if (length(condition_labels) != 2) stop("condition_labels must have exactly two elements: c(control, stress)")
  if (is.null(title)) {
    title <- sprintf("Expression models (%s = control, %s = stress)", condition_labels[1], condition_labels[2])
    if (!is.null(classification)) title <- paste0(title, "; % of responsive features")
  }
  tick_cex <- if (max(nchar(condition_labels)) <= 8) 0.7 else 0.6
  if (is.null(group_names)) {
    group_names <- if (!is.null(classification) && !is.null(attr(classification, "group_names"))) {
      attr(classification, "group_names")
    } else c("A", "B")
  }
  summ <- if (!is.null(classification)) summarize_classification(classification) else NULL
  op <- graphics::par(mfrow = c(5, 4), mar = c(2.2, 1, 1.6, 0.5), oma = c(2, 3, 2, 0))
  on.exit(graphics::par(op))
  for (model in names(.patterns)) {
    pats <- .patterns[[model]]
    cols <- .columns[[model]]
    for (col in 1:4) {
      j <- match(col, cols)
      if (is.na(j)) {
        graphics::plot.new()
        next
      }
      lv <- pats[[j]]
      graphics::plot(NA, xlim = c(0, 6), ylim = c(-1.3, 1.3), axes = FALSE, xlab = "", ylab = "")
      graphics::box(col = "grey60")
      graphics::abline(v = c(1.5, 4.5), lty = 2, col = "grey70")
      graphics::abline(v = 3, col = "grey70")
      graphics::abline(h = 0, lty = 2, col = "grey70")
      graphics::lines(c(0, 1, 2, 3), lv[c(1, 1, 2, 2)], col = .col_a, lwd = 3)
      graphics::lines(c(3, 4, 5, 6), lv[c(3, 3, 4, 4)], col = .col_b, lwd = 3)
      graphics::axis(3, at = c(0.75, 2.25, 3.75, 5.25), labels = rep(condition_labels, 2),
                     tick = FALSE, line = -0.9, cex.axis = tick_cex)
      lab <- paste("Case", j)
      if (!is.null(summ)) {
        r <- summ[summ$model == model & summ$case == j, ]
        n <- if (nrow(r)) r$n else 0
        pct <- if (nrow(r)) r$pct_of_responsive else NA
        lab <- paste0(lab, ": ", n, if (!is.na(pct)) sprintf(" (%.1f%%)", pct) else "")
      }
      graphics::mtext(lab, side = 1, line = 0.6, cex = 0.6)
      if (col == 1) graphics::mtext(row_labels[[model]], side = 2, line = 0.8, cex = 0.65)
    }
  }
  graphics::mtext(title, outer = TRUE, side = 3, cex = 0.8)
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot.new()
  graphics::legend("bottom", legend = group_names, col = c(.col_a, .col_b), lwd = 3, horiz = TRUE, bty = "n",
                   cex = 0.8)
  invisible(NULL)
}

#' Bar chart of features per model and case
#'
#' @param classification Output of [classify_responses()].
#' @param ... Passed to [graphics::barplot()].
#' @return Invisibly, the summary table.
#' @export
plot_summary <- function(classification, ...) {
  s <- summarize_classification(classification)
  op <- graphics::par(mar = c(4, 9, 1, 2))
  on.exit(graphics::par(op))
  bp <- graphics::barplot(rev(s$n), names.arg = rev(s$label), horiz = TRUE, las = 1, col = "#4c72b0",
                          xlab = "Features", cex.names = 0.75, ...)
  graphics::text(rev(s$n), bp, labels = rev(s$n), pos = 4, cex = 0.7, xpd = TRUE)
  invisible(s)
}
