.patterns <- list(
  stress_ready = list(c(0, 0, -1, 0), c(0, 0, 1, 0), c(0, 1, 1, 1), c(0, -1, -1, -1)),
  shared = list(c(0, 1, 0, 1), c(0, -1, 0, -1)),
  unique = list(c(0, 1, 0, 0), c(0, -1, 0, 0), c(0, 0, 0, 1), c(0, 0, 0, -1)),
  opposite = list(c(0, 1, 0, -1), c(0, -1, 0, 1)),
  no_response = list(c(0, 0, 1, 1), c(0, 0, -1, -1), c(0, 0, 0, 0))
)
.titles <- c(stress_ready = "Stress-ready", shared = "Shared response", unique = "Unique response",
             opposite = "Opposite response", no_response = "No response")
.col_a <- "#1f77b4"
.col_b <- "#ff7f0e"

#' Plot the expression model patterns
#'
#' Draws the 15 schematic patterns (C = control, S = stress; left half group A,
#' right half group B), annotated with counts and percentages of responsive
#' features when a classification is supplied. Uses base graphics.
#'
#' @param classification Optional output of [classify_responses()].
#' @param group_names Names for the legend; taken from `classification` when
#'   available.
#' @return Invisibly, `NULL`.
#' @examples
#' plot_model_patterns()
#' @export
plot_model_patterns <- function(classification = NULL, group_names = NULL) {
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
    for (j in 1:4) {
      if (j > length(pats)) {
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
      graphics::axis(3, at = c(0.75, 2.25, 3.75, 5.25), labels = c("C", "S", "C", "S"),
                     tick = FALSE, line = -0.9, cex.axis = 0.7)
      lab <- paste("Case", j)
      if (!is.null(summ)) {
        r <- summ[summ$model == model & summ$case == j, ]
        n <- if (nrow(r)) r$n else 0
        pct <- if (nrow(r)) r$pct_of_responsive else NA
        lab <- paste0(lab, ": ", n, if (!is.na(pct)) sprintf(" (%.1f%%)", pct) else "")
      }
      graphics::mtext(lab, side = 1, line = 0.6, cex = 0.6)
      if (j == 1) graphics::mtext(.titles[[model]], side = 2, line = 0.8, cex = 0.65)
    }
  }
  graphics::mtext("Expression models (C = control, S = stress)", outer = TRUE, side = 3, cex = 0.8)
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
