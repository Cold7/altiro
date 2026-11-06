#' Contrast keys used by the classifier
#'
#' The classification rules need five differential-expression contrasts between the
#' four cells of a 2 x 2 design (groups A and B, conditions control and stress):
#'
#' \describe{
#'   \item{a_response}{A stress vs A control}
#'   \item{b_response}{B stress vs B control}
#'   \item{basal}{B control vs A control}
#'   \item{b_stress_vs_a_control}{B stress vs A control}
#'   \item{a_stress_vs_b_control}{A stress vs B control}
#' }
#'
#' @return A named character vector (keys as names, descriptions as values).
#' @examples
#' contrast_keys()
#' @export
contrast_keys <- function() {
  c(
    a_response = "A stress vs A control",
    b_response = "B stress vs B control",
    basal = "B control vs A control",
    b_stress_vs_a_control = "B stress vs A control",
    a_stress_vs_b_control = "A stress vs B control"
  )
}

# "stress_ready" is the legacy name of "primed" (label_style = "stress_ready")
.primed_models <- c("primed", "stress_ready")
.responsive_models <- c("primed", "stress_ready", "shared", "unique", "opposite", "unresolved")
.model_order <- c("primed", "stress_ready", "shared", "unique", "opposite", "no_response", "unresolved", "unclassified")

.definitions <- function() {
  data.frame(
    model = c(rep("primed", 4), rep("shared", 2), rep("unique", 4), rep("opposite", 2),
              rep("no_response", 3), rep("unresolved", 4), "unclassified"),
    case = c(1:4, 1:2, 1:4, 1:2, 1:3, 1:4, 0L),
    primed = c("A", "A", "B", "B", rep("", 16)),
    responding = c("B", "B", "A", "A", "both", "both", "A", "A", "B", "B", "both", "both",
                   "none", "none", "none", "B", "B", "A", "A", ""),
    text = c(
      "{B} is induced under stress up to the control level of {A} ({A} is primed)",
      "{B} is repressed under stress down to the control level of {A} ({A} is primed)",
      "{A} is induced under stress up to the control level of {B} ({B} is primed)",
      "{A} is repressed under stress down to the control level of {B} ({B} is primed)",
      "Induced under stress in both {A} and {B}",
      "Repressed under stress in both {A} and {B}",
      "Induced under stress only in {A}",
      "Repressed under stress only in {A}",
      "Induced under stress only in {B}",
      "Repressed under stress only in {B}",
      "Induced in {A} and repressed in {B}",
      "Repressed in {A} and induced in {B}",
      "No response; basal level higher in {B} than in {A}",
      "No response; basal level lower in {B} than in {A}",
      "No response; no basal difference detected",
      rep("Responds only in {B}; level under stress neither different from nor equivalent to {A} control", 2),
      rep("Responds only in {A}; level under stress neither different from nor equivalent to {B} control", 2),
      "Missing statistics in a contrast required by the classification rules"
    ),
    stringsAsFactors = FALSE
  )
}

#' Definitions of the expression models
#'
#' Table describing every model and case produced by [classify_responses()].
#'
#' @param group_names Character vector of length 2 with display names for
#'   groups A and B.
#' @param include_extra Also list `unresolved` and `unclassified`, which are
#'   only produced by non-default options.
#' @return A data frame with columns `model`, `case`, `label`, `primed_group`,
#'   `responding_group` and `description`.
#' @examples
#' model_definitions(c("S. lycopersicum", "S. pennellii"))
#' @export
model_definitions <- function(group_names = c("A", "B"), include_extra = FALSE) {
  stopifnot(length(group_names) == 2)
  d <- .definitions()
  if (!include_extra) d <- d[!d$model %in% c("unresolved", "unclassified"), ]
  map <- function(x) ifelse(x == "A", group_names[1], ifelse(x == "B", group_names[2], x))
  desc <- gsub("{A}", group_names[1], d$text, fixed = TRUE)
  desc <- gsub("{B}", group_names[2], desc, fixed = TRUE)
  out <- data.frame(
    model = d$model,
    case = as.integer(d$case),
    label = ifelse(d$model == "unclassified", "unclassified", paste0(d$model, "_", d$case)),
    primed_group = map(d$primed),
    responding_group = map(d$responding),
    description = desc,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

#' Classification options and presets
#'
#' @param preset `"paper"` (default) reproduces exactly the classification rules of
#'   Contreras-Riquelme et al. (2025). `"strict"` sets
#'   `na_action = "unclassified"` and `require_basal_difference = TRUE`.
#' @param alpha Adjusted p-value threshold (default 0.05).
#' @param min_abs_lfc Minimum absolute log2 fold change required, in addition
#'   to `alpha`, to call a contrast significant (default 0).
#' @param na_action `"not_de"`: a missing adjusted p-value counts as not
#'   differentially expressed (paper behaviour). `"unclassified"`: features
#'   with a missing value in any contrast consulted by the tree are labelled
#'   `unclassified`.
#' @param require_basal_difference If `TRUE`, primed calls also require
#'   a significant basal difference in the direction implied by the case;
#'   otherwise the feature is labelled as the corresponding unique response.
#' @param equivalence_margin `NULL` (default) or a positive number (log2
#'   units). If set, "same level" must be shown with a two one-sided test
#'   ([equivalence_test()]) on the matching contrast; features neither
#'   different nor equivalent are labelled `unresolved`. Needs `lfcSE`.
#' @param sign_digits Decimals used to round log2 fold changes before taking
#'   their sign (the published script used 2).
#' @return A list of class `altiro_options`.
#' @examples
#' classification_options("strict")
#' classification_options("paper", equivalence_margin = 0.5)
#' @export
classification_options <- function(preset = c("paper", "strict"), alpha = NULL, min_abs_lfc = NULL,
                                   na_action = NULL, require_basal_difference = NULL,
                                   equivalence_margin = NULL, sign_digits = NULL) {
  preset <- match.arg(preset)
  o <- list(alpha = 0.05, min_abs_lfc = 0, na_action = "not_de", require_basal_difference = FALSE,
            equivalence_margin = NULL, sign_digits = 2L)
  if (preset == "strict") {
    o$na_action <- "unclassified"
    o$require_basal_difference <- TRUE
  }
  if (!is.null(alpha)) o$alpha <- alpha
  if (!is.null(min_abs_lfc)) o$min_abs_lfc <- min_abs_lfc
  if (!is.null(na_action)) o$na_action <- na_action
  if (!is.null(require_basal_difference)) o$require_basal_difference <- require_basal_difference
  if (!is.null(equivalence_margin)) o$equivalence_margin <- equivalence_margin
  if (!is.null(sign_digits)) o$sign_digits <- as.integer(sign_digits)
  if (!(o$alpha > 0 && o$alpha < 1)) stop("alpha must be in (0, 1)")
  if (o$min_abs_lfc < 0) stop("min_abs_lfc must be >= 0")
  if (!o$na_action %in% c("not_de", "unclassified")) stop("na_action must be 'not_de' or 'unclassified'")
  if (!is.null(o$equivalence_margin) && o$equivalence_margin <= 0) stop("equivalence_margin must be > 0")
  o$preset <- preset
  class(o) <- "altiro_options"
  o
}

#' @export
print.altiro_options <- function(x, ...) {
  cat("altiro classification options (preset:", x$preset, ")\n")
  cat("  alpha                    :", x$alpha, "\n")
  cat("  min_abs_lfc              :", x$min_abs_lfc, "\n")
  cat("  na_action                :", x$na_action, "\n")
  cat("  require_basal_difference :", x$require_basal_difference, "\n")
  cat("  equivalence_margin       :", if (is.null(x$equivalence_margin)) "none" else x$equivalence_margin, "\n")
  cat("  sign_digits              :", x$sign_digits, "\n")
  invisible(x)
}

#' Equivalence test for log2 fold changes (TOST)
#'
#' Two one-sided Wald tests of `|log2FC| < margin`. The p-value is the larger
#' of the two one-sided p-values, as in DESeq2 with `altHypothesis = "lessAbs"`.
#'
#' @param log2fc Numeric vector of log2 fold changes.
#' @param lfc_se Numeric vector of their standard errors.
#' @param margin Equivalence margin in log2 units.
#' @return A data frame with `pvalue` and Benjamini-Hochberg `padj`.
#' @examples
#' equivalence_test(c(0, 0.3, 1.2), c(0.1, 0.2, 0.2), margin = 0.5)
#' @export
equivalence_test <- function(log2fc, lfc_se, margin) {
  p_low <- stats::pnorm((log2fc + margin) / lfc_se, lower.tail = FALSE)
  p_high <- stats::pnorm((log2fc - margin) / lfc_se)
  p <- pmax(p_low, p_high)
  p[is.na(log2fc) | is.na(lfc_se) | lfc_se <= 0] <- NA_real_
  data.frame(pvalue = p, padj = stats::p.adjust(p, method = "BH"), row.names = names(log2fc))
}

.as_tables <- function(contrasts) {
  if (inherits(contrasts, "contrast_set")) contrasts <- contrasts$tables
  keys <- names(contrast_keys())
  if (is.data.frame(contrasts)) {
    tabs <- list()
    for (k in keys) {
      cols <- grep(paste0("^", k, "_"), colnames(contrasts), value = TRUE)
      if (!length(cols)) stop("No columns starting with '", k, "_' in the wide table")
      d <- contrasts[, cols, drop = FALSE]
      colnames(d) <- substring(cols, nchar(k) + 2)
      tabs[[k]] <- d
    }
    contrasts <- tabs
  }
  missing <- setdiff(keys, names(contrasts))
  if (length(missing)) stop("Missing contrasts: ", paste(missing, collapse = ", "))
  out <- list()
  for (k in keys) {
    d <- as.data.frame(contrasts[[k]])
    for (col in c("log2FoldChange", "padj")) {
      if (!col %in% colnames(d)) stop("Contrast '", k, "' lacks column '", col, "'")
    }
    if (anyDuplicated(rownames(d))) stop("Contrast '", k, "' has duplicated feature identifiers")
    out[[k]] <- d
  }
  out
}

.sign_rounded <- function(lfc, digits) {
  # Same rule as numpy.round(x, digits): round half to even on x * 10^digits.
  r <- round(lfc * 10^digits)
  s <- ifelse(r > 0, 1L, ifelse(r < 0, -1L, 0L))
  s[is.na(s)] <- 0L
  s
}

#' Classify features into expression response models
#'
#' Applies the classification rules of Contreras-Riquelme et al. (2025) to the five
#' contrasts of a 2 x 2 design and assigns each feature (ortholog pair or gene)
#' to one model: primed (cases 1-4; called stress-ready in the original
#' publication), shared (1-2), unique (1-4),
#' opposite (1-2) or no response (1-3). With the default `preset = "paper"`
#' the result is identical to the published script, where A =
#' *S. lycopersicum* and B = *S. pennellii*.
#'
#' @param contrasts A `contrast_set` from [run_contrasts()], a named list with
#'   the five data frames listed in [contrast_keys()] (feature ids as row
#'   names; columns `log2FoldChange`, `padj` and, for equivalence tests,
#'   `lfcSE`), or one wide data frame with columns `<key>_log2FoldChange`,
#'   `<key>_padj`, ...
#' @param preset `"paper"` or `"strict"`; see [classification_options()].
#' @param group_names Display names for groups A and B.
#' @param options An object from [classification_options()]; overrides
#'   `preset` and `...`.
#' @param keep_statistics Append `<key>_log2FoldChange` and `<key>_padj`.
#' @param label_style `"primed"` (default) or `"stress_ready"`, the name used
#'   before version 0.5.0 and in the original publication. Only the names in
#'   the `label` and `model` columns change; the classification is the same.
#' @param ... Individual options passed to [classification_options()], e.g.
#'   `alpha = 0.01` or `equivalence_margin = 0.5`.
#' @return A data frame (one row per feature, row names = feature ids) with
#'   columns `feature`, `label`, `model`, `case`, `primed_group`,
#'   `responding_group`, `description`, `flags` and, optionally, the
#'   statistics of each contrast.
#' @examples
#' f <- function(k) system.file("extdata", paste0("contrast_", k, ".tsv"), package = "altiro")
#' paths <- setNames(lapply(names(contrast_keys()), f), names(contrast_keys()))
#' tabs <- read_contrasts(paths)
#' res <- classify_responses(tabs, group_names = c("sensitive", "tolerant"))
#' table(res$label)
#' @export
classify_responses <- function(contrasts, preset = c("paper", "strict"), group_names = c("A", "B"),
                               options = NULL, keep_statistics = TRUE,
                               label_style = c("primed", "stress_ready"), ...) {
  preset <- match.arg(preset)
  label_style <- match.arg(label_style)
  opts <- if (is.null(options)) classification_options(preset, ...) else options
  tabs <- .as_tables(contrasts)
  keys <- names(contrast_keys())
  ids <- rownames(tabs[[1]])
  sizes <- vapply(tabs, nrow, integer(1))
  for (k in keys[-1]) ids <- ids[ids %in% rownames(tabs[[k]])]
  if (length(unique(sizes)) > 1 || length(ids) != sizes[1]) {
    warning("Contrast tables cover different features; using the ", length(ids), " shared ones.")
  }
  tabs <- lapply(tabs, function(d) d[ids, , drop = FALSE])
  n <- length(ids)

  st <- function(k) {
    padj <- suppressWarnings(as.numeric(tabs[[k]]$padj))
    lfc <- suppressWarnings(as.numeric(tabs[[k]]$log2FoldChange))
    na <- is.na(padj)
    lfc0 <- ifelse(is.na(lfc), 0, lfc)
    de <- !na & padj < opts$alpha & abs(lfc0) >= opts$min_abs_lfc
    de[is.na(de)] <- FALSE
    list(de = de, na = na, s = .sign_rounded(lfc, opts$sign_digits), lfc = lfc, padj = padj)
  }
  A <- st("a_response"); B <- st("b_response"); BAS <- st("basal")
  MB <- st("b_stress_vs_a_control"); MA <- st("a_stress_vs_b_control")

  if (!is.null(opts$equivalence_margin)) {
    for (k in c("b_stress_vs_a_control", "a_stress_vs_b_control")) {
      if (!"lfcSE" %in% colnames(tabs[[k]])) stop("equivalence_margin requires an 'lfcSE' column")
    }
    eqb <- equivalence_test(MB$lfc, as.numeric(tabs$b_stress_vs_a_control$lfcSE), opts$equivalence_margin)$padj
    eqa <- equivalence_test(MA$lfc, as.numeric(tabs$a_stress_vs_b_control$lfcSE), opts$equivalence_margin)$padj
    same_b <- !MB$de & !is.na(eqb) & eqb < opts$alpha
    same_a <- !MA$de & !is.na(eqa) & eqa < opts$alpha
  } else {
    same_b <- !MB$de
    same_a <- !MA$de
  }

  model <- rep("", n)
  case <- integer(n)
  flags <- vector("list", n)

  only_b <- !A$de & B$de
  only_a <- A$de & !B$de
  both <- A$de & B$de
  none <- !A$de & !B$de
  b_up <- B$s == 1L
  a_up <- A$s == 1L

  sel <- only_b & same_b; model[sel] <- "primed"; case[sel] <- ifelse(b_up, 1L, 2L)[sel]
  sel <- only_b & MB$de; model[sel] <- "unique"; case[sel] <- ifelse(b_up, 3L, 4L)[sel]
  sel <- only_b & !MB$de & !same_b; model[sel] <- "unresolved"; case[sel] <- ifelse(b_up, 1L, 2L)[sel]

  sel <- only_a & same_a; model[sel] <- "primed"; case[sel] <- ifelse(a_up, 3L, 4L)[sel]
  sel <- only_a & MA$de; model[sel] <- "unique"; case[sel] <- ifelse(a_up, 1L, 2L)[sel]
  sel <- only_a & !MA$de & !same_a; model[sel] <- "unresolved"; case[sel] <- ifelse(a_up, 3L, 4L)[sel]

  sel <- both & A$s == 1L & B$s == 1L; model[sel] <- "shared"; case[sel] <- 1L
  sel <- both & A$s == -1L & B$s == -1L; model[sel] <- "shared"; case[sel] <- 2L
  sel <- both & A$s == 1L & B$s == -1L; model[sel] <- "opposite"; case[sel] <- 1L
  sel <- both & model == ""; model[sel] <- "opposite"; case[sel] <- 2L

  sel <- none & BAS$de & BAS$s == 1L; model[sel] <- "no_response"; case[sel] <- 1L
  sel <- none & BAS$de & BAS$s == -1L; model[sel] <- "no_response"; case[sel] <- 2L
  sel <- none & model == ""; model[sel] <- "no_response"; case[sel] <- 3L

  add_flag <- function(idx, flag) for (i in which(idx)) flags[[i]] <<- c(flags[[i]], flag)

  if (isTRUE(opts$require_basal_difference)) {
    sr <- model == "primed"
    expected <- c(-1L, 1L, 1L, -1L)[pmax(case, 1L)]
    ok <- BAS$de & BAS$s == expected
    demote <- sr & !ok
    case[demote] <- c(3L, 4L, 1L, 2L)[case[demote]]
    model[demote] <- "unique"
    add_flag(demote, "basal_difference_not_significant")
  }

  consulted_na <- A$na | B$na | (only_b & MB$na) | (only_a & MA$na) | (none & BAS$na)
  if (isTRUE(opts$require_basal_difference)) consulted_na <- consulted_na | ((only_a | only_b) & BAS$na)
  if (opts$na_action == "unclassified") {
    model[consulted_na] <- "unclassified"
    case[consulted_na] <- 0L
  } else {
    add_flag(consulted_na, "na_treated_as_not_de")
  }

  defs <- model_definitions(group_names, include_extra = TRUE)
  key <- paste(model, case)
  m <- match(key, paste(defs$model, defs$case))
  if (label_style != "primed") model[model == "primed"] <- label_style
  out <- data.frame(
    feature = ids,
    label = ifelse(model == "unclassified", "unclassified", paste0(model, "_", case)),
    model = model,
    case = case,
    primed_group = defs$primed_group[m],
    responding_group = defs$responding_group[m],
    description = defs$description[m],
    flags = vapply(flags, function(f) paste(f, collapse = ";"), character(1)),
    stringsAsFactors = FALSE
  )
  if (keep_statistics) {
    for (k in keys) {
      out[[paste0(k, "_log2FoldChange")]] <- suppressWarnings(as.numeric(tabs[[k]]$log2FoldChange))
      out[[paste0(k, "_padj")]] <- suppressWarnings(as.numeric(tabs[[k]]$padj))
    }
  }
  rownames(out) <- ids
  attr(out, "options") <- opts
  attr(out, "group_names") <- group_names
  attr(out, "label_style") <- label_style
  out
}

#' Summarise a classification
#'
#' Counts per model and case. `pct_of_responsive` uses as denominator the
#' features that respond in at least one group (primed, shared, unique,
#' opposite and unresolved); `pct_of_all` uses every feature.
#'
#' @param classification Output of [classify_responses()].
#' @return A data frame with `label`, `model`, `case`, `n`,
#'   `pct_of_responsive` and `pct_of_all`.
#' @export
summarize_classification <- function(classification) {
  tab <- stats::aggregate(list(n = rep(1L, nrow(classification))),
                          by = list(model = classification$model, case = classification$case), FUN = sum)
  tab <- tab[order(match(tab$model, .model_order), tab$case), ]
  tab$label <- ifelse(tab$model == "unclassified", "unclassified", paste0(tab$model, "_", tab$case))
  n_resp <- sum(classification$model %in% .responsive_models)
  tab$pct_of_responsive <- ifelse(tab$model %in% .responsive_models, 100 * tab$n / max(n_resp, 1), NA_real_)
  tab$pct_of_all <- 100 * tab$n / max(nrow(classification), 1)
  tab <- tab[, c("label", "model", "case", "n", "pct_of_responsive", "pct_of_all")]
  rownames(tab) <- NULL
  tab
}

#' Test whether primed states are more frequent in one direction
#'
#' Among features that respond to stress in only one group, compares the
#' fraction that reaches the other group's control level (primed) for
#' A-only and B-only responders with Fisher's exact test. A larger raw count of
#' primed genes in one direction can simply reflect more responsive genes
#' in that group; this comparison conditions on it.
#'
#' @param classification Output of [classify_responses()].
#' @return A data frame with one row per responding group, the counts,
#'   `fraction_primed`, the sample `odds_ratio` and `fisher_pvalue`.
#' @export
primed_asymmetry <- function(classification) {
  m <- classification$model
  cs <- classification$case
  a_sr <- sum(m %in% .primed_models & cs %in% 3:4)
  a_un <- sum(m == "unique" & cs %in% 1:2)
  b_sr <- sum(m %in% .primed_models & cs %in% 1:2)
  b_un <- sum(m == "unique" & cs %in% 3:4)
  gn <- attr(classification, "group_names")
  if (is.null(gn)) gn <- c("A", "B")
  mat <- matrix(c(a_sr, b_sr, a_un, b_un), nrow = 2)
  p <- stats::fisher.test(mat)$p.value
  or <- if (a_un * b_sr == 0) Inf else (a_sr * b_un) / (a_un * b_sr)
  frac <- function(sr, un) if (sr + un == 0) NA_real_ else sr / (sr + un)
  data.frame(
    responding_group = gn, primed_group = rev(gn),
    primed = c(a_sr, b_sr), unique = c(a_un, b_un),
    fraction_primed = c(frac(a_sr, a_un), frac(b_sr, b_un)),
    odds_ratio = or, fisher_pvalue = p,
    stringsAsFactors = FALSE
  )
}
