.cells <- c("A_control", "A_stress", "B_control", "B_stress")
.contrast_cells <- list(
  a_response = c("A_stress", "A_control"),
  b_response = c("B_stress", "B_control"),
  basal = c("B_control", "A_control"),
  b_stress_vs_a_control = c("B_stress", "A_control"),
  a_stress_vs_b_control = c("A_stress", "B_control")
)

#' Map sample annotations to the four design cells
#'
#' @param samples Data frame of sample annotations; row names are sample names.
#' @param group_col,condition_col Columns holding the group and the condition.
#' @param group_a,group_b Levels of `group_col`: A is the reference (e.g. the
#'   sensitive species or line), B the comparison (e.g. the tolerant one).
#' @param control,stress Levels of `condition_col`.
#' @return A data frame with columns `group`, `condition` and `cell`
#'   (`A_control`, `A_stress`, `B_control`, `B_stress`). Samples outside the
#'   design are dropped with a warning.
#' @export
build_design <- function(samples, group_col, condition_col, group_a, group_b, control, stress) {
  for (col in c(group_col, condition_col)) {
    if (!col %in% colnames(samples)) stop("Column '", col, "' not found in the sample table")
  }
  g <- as.character(samples[[group_col]])
  cnd <- as.character(samples[[condition_col]])
  keep <- g %in% c(group_a, group_b) & cnd %in% c(control, stress)
  if (any(!keep)) {
    warning("Dropping ", sum(!keep), " sample(s) outside the 2 x 2 design: ",
            paste(rownames(samples)[!keep], collapse = ", "))
  }
  d <- data.frame(group = g[keep], condition = cnd[keep], row.names = rownames(samples)[keep],
                  stringsAsFactors = FALSE)
  d$cell <- paste0(ifelse(d$group == group_a, "A", "B"), "_", ifelse(d$condition == control, "control", "stress"))
  sizes <- table(factor(d$cell, levels = .cells))
  if (any(sizes == 0)) stop("Empty design cell(s): ", paste(names(sizes)[sizes == 0], collapse = ", "))
  if (any(sizes < 2)) warning("Some design cells have fewer than 2 replicates")
  d$cell <- factor(d$cell, levels = .cells)
  d
}

#' Run DESeq2 and extract the five contrasts
#'
#' Fits one DESeq2 model with a four-level factor (`~ cell`) over all samples,
#' as in the original publication, and extracts the contrasts listed in
#' [contrast_keys()].
#'
#' @param counts Raw count matrix or data frame (features in rows, samples in
#'   columns).
#' @param samples Sample annotation; row names must match `colnames(counts)`.
#' @inheritParams build_design
#' @param min_total_count Features with a total count `<=` this value are
#'   removed before fitting (the published script used `rowSums > 10`).
#' @param alpha Significance level used by the independent filtering unless
#'   `filter_alpha` is given. The published script filtered with 0.01 and
#'   classified with 0.05.
#' @param filter_alpha Optional separate level for independent filtering.
#' @param cooks_filter,independent_filter Passed to [DESeq2::results()].
#' @param round_counts Round non-integer input (e.g. estimated counts).
#' @param quiet Silence DESeq2 messages.
#' @return An object of class `contrast_set`: a list with `tables` (named list
#'   of five data frames), `normalized_counts`, `design`, `group_names` and
#'   `info`.
#' @export
run_contrasts <- function(counts, samples, group_col, condition_col, group_a, group_b, control, stress,
                          min_total_count = 10, alpha = 0.05, filter_alpha = NULL, cooks_filter = TRUE,
                          independent_filter = TRUE, round_counts = FALSE, quiet = TRUE) {
  design <- build_design(samples, group_col, condition_col, group_a, group_b, control, stress)
  missing <- setdiff(rownames(design), colnames(counts))
  if (length(missing)) stop("Samples missing from the count matrix: ", paste(missing, collapse = ", "))
  mat <- as.matrix(counts[, rownames(design), drop = FALSE])
  storage.mode(mat) <- "double"
  if (anyNA(mat)) stop("The count matrix contains missing values")
  if (any(abs(mat - round(mat)) > 1e-8)) {
    if (!round_counts) stop("Counts are not integers. Use raw counts, or set round_counts = TRUE.")
    mat <- round(mat)
  }
  if (any(mat < 0)) stop("Counts must be non-negative")
  if (anyDuplicated(rownames(mat))) stop("Duplicated feature identifiers in the count matrix")
  storage.mode(mat) <- "integer"
  n_input <- nrow(mat)
  mat <- mat[rowSums(mat) > min_total_count, , drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = design, design = ~cell)
  dds <- DESeq2::DESeq(dds, quiet = quiet)
  f_alpha <- if (is.null(filter_alpha)) alpha else filter_alpha
  tables <- list()
  for (k in names(contrast_keys())) {
    cc <- .contrast_cells[[k]]
    args <- list(object = dds, contrast = c("cell", cc[1], cc[2]), alpha = f_alpha,
                 independentFiltering = independent_filter, pAdjustMethod = "BH")
    if (!cooks_filter) args$cooksCutoff <- FALSE
    res <- do.call(DESeq2::results, args)
    tables[[k]] <- as.data.frame(res)
  }
  norm <- BiocGenerics::counts(dds, normalized = TRUE)
  structure(
    list(
      tables = tables,
      normalized_counts = as.data.frame(norm),
      design = design,
      group_names = c(as.character(group_a), as.character(group_b)),
      info = list(n_input_features = n_input, n_tested_features = nrow(mat),
                  min_total_count = min_total_count, filter_alpha = f_alpha, engine = "DESeq2")
    ),
    class = "contrast_set"
  )
}

#' @export
print.contrast_set <- function(x, ...) {
  cat("contrast_set (", x$info$engine, "): ", x$info$n_tested_features, " features tested; groups A = ",
      x$group_names[1], ", B = ", x$group_names[2], "\n", sep = "")
  for (k in names(x$tables)) {
    cat(sprintf("  %-22s %-24s significant (padj < 0.05): %d\n", k, contrast_keys()[[k]],
                sum(x$tables[[k]]$padj < 0.05, na.rm = TRUE)))
  }
  invisible(x)
}

#' Read five contrast tables from files
#'
#' @param paths Named list or vector (names from [contrast_keys()]) of file
#'   paths. For the files written by the original publication script use
#'   `list(a_response = "Sly.drought_Sly.control.csv",
#'   b_response = "Spenn.drought_Spenn.control.csv",
#'   basal = "Spenn.control_Sly.control.csv",
#'   b_stress_vs_a_control = "Spenn.drought_Sly.control.csv",
#'   a_stress_vs_b_control = "Sly.drought_Spenn.control.csv")` with
#'   `id_col = "Row.names"`.
#' @param id_col Column holding the feature identifier (default: first column).
#' @param sep Field separator; inferred from the extension when `NULL`.
#' @param lfc_col,padj_col,se_col Column names mapped onto `log2FoldChange`,
#'   `padj` and `lfcSE`.
#' @return A named list of five data frames with feature ids as row names.
#' @export
read_contrasts <- function(paths, id_col = NULL, sep = NULL, lfc_col = "log2FoldChange",
                           padj_col = "padj", se_col = "lfcSE") {
  out <- list()
  for (k in names(contrast_keys())) {
    if (is.null(paths[[k]])) stop("No path given for contrast '", k, "'")
    p <- paths[[k]]
    s <- if (!is.null(sep)) sep else if (grepl("\\.(tsv|txt|tab)(\\.gz)?$", p, ignore.case = TRUE)) "\t" else ","
    d <- utils::read.delim(p, sep = s, check.names = FALSE, stringsAsFactors = FALSE)
    col <- if (is.null(id_col)) colnames(d)[1] else id_col
    if (!col %in% colnames(d)) stop(p, ": id column '", col, "' not found")
    rownames(d) <- as.character(d[[col]])
    d[[col]] <- NULL
    names(d)[names(d) == lfc_col] <- "log2FoldChange"
    names(d)[names(d) == padj_col] <- "padj"
    names(d)[names(d) == se_col] <- "lfcSE"
    out[[k]] <- d
  }
  out
}

#' Write contrast tables to a directory
#'
#' @param contrasts A `contrast_set` or named list of five data frames.
#' @param outdir Output directory; files are named `contrast_<key>.tsv`.
#' @return Invisibly, the written paths.
#' @export
write_contrasts <- function(contrasts, outdir) {
  tabs <- if (inherits(contrasts, "contrast_set")) contrasts$tables else contrasts
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  paths <- character()
  for (k in names(contrast_keys())) {
    p <- file.path(outdir, paste0("contrast_", k, ".tsv"))
    d <- data.frame(feature = rownames(tabs[[k]]), tabs[[k]], check.names = FALSE)
    utils::write.table(d, p, sep = "\t", quote = FALSE, row.names = FALSE)
    paths[k] <- p
  }
  invisible(paths)
}
