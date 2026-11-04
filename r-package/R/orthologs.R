.blast6 <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
             "qstart", "qend", "sstart", "send", "evalue", "bitscore")

#' Read BLAST tabular output
#'
#' Reads BLAST `-outfmt 6` output (12 standard columns).
#'
#' @param path File path.
#' @param strip_pattern Optional regular expression removed from query and
#'   subject identifiers, e.g. `"\\.\\d+$"` to drop a version suffix.
#' @return A data frame with the 12 standard column names.
#' @export
read_blast_tabular <- function(path, strip_pattern = NULL) {
  d <- utils::read.delim(path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
  if (ncol(d) < 12) stop(path, ": expected at least 12 columns (BLAST -outfmt 6)")
  d <- d[, 1:12]
  colnames(d) <- .blast6
  .strip(d, strip_pattern)
}

.strip <- function(d, strip_pattern) {
  d$qseqid <- as.character(d$qseqid)
  d$sseqid <- as.character(d$sseqid)
  if (!is.null(strip_pattern)) {
    d$qseqid <- sub(strip_pattern, "", d$qseqid, perl = TRUE)
    d$sseqid <- sub(strip_pattern, "", d$sseqid, perl = TRUE)
  }
  d
}

#' Best hit per query
#'
#' Highest bitscore, then lowest e-value, then subject identifier (deterministic
#' tie-break).
#'
#' @param blast Data frame from [read_blast_tabular()].
#' @return One row per query.
#' @export
best_hits <- function(blast) {
  o <- order(blast$qseqid, -blast$bitscore, blast$evalue, blast$sseqid, method = "radix")
  b <- blast[o, ]
  b <- b[!duplicated(b$qseqid), ]
  rownames(b) <- NULL
  b
}

#' Orthologs as reciprocal best BLAST hits
#'
#' @param blast_ab BLAST of proteome A (queries) against proteome B: a path or
#'   a data frame with the standard column names.
#' @param blast_ba BLAST of proteome B against proteome A.
#' @param max_evalue,min_bitscore Optional filters applied before choosing
#'   best hits.
#' @param strip_pattern Regular expression removed from identifiers.
#' @return A data frame with `ortholog_id`, `gene_a`, `gene_b`, `pident_ab`,
#'   `evalue_ab`, `bitscore_ab`, `evalue_ba` and `bitscore_ba`.
#' @export
reciprocal_best_hits <- function(blast_ab, blast_ba, max_evalue = NULL, min_bitscore = NULL,
                                 strip_pattern = NULL) {
  load <- function(x) if (is.data.frame(x)) .strip(x, strip_pattern) else read_blast_tabular(x, strip_pattern)
  ab <- load(blast_ab)
  ba <- load(blast_ba)
  filt <- function(d) {
    if (!is.null(max_evalue)) d <- d[d$evalue <= max_evalue, ]
    if (!is.null(min_bitscore)) d <- d[d$bitscore >= min_bitscore, ]
    d
  }
  bab <- best_hits(filt(ab))
  bba <- best_hits(filt(ba))
  key_ab <- paste(bab$qseqid, bab$sseqid, sep = "\r")
  key_ba <- paste(bba$sseqid, bba$qseqid, sep = "\r")
  i <- match(key_ab, key_ba)
  keep <- !is.na(i)
  out <- data.frame(
    gene_a = bab$qseqid[keep], gene_b = bab$sseqid[keep], pident_ab = bab$pident[keep],
    evalue_ab = bab$evalue[keep], bitscore_ab = bab$bitscore[keep],
    evalue_ba = bba$evalue[i[keep]], bitscore_ba = bba$bitscore[i[keep]],
    stringsAsFactors = FALSE
  )
  out <- out[order(out$gene_a, method = "radix"), ]
  out <- data.frame(ortholog_id = paste0("ortholog_", seq_len(nrow(out))), out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' Merge two species' count matrices through a one-to-one ortholog table
#'
#' @param counts_a,counts_b Raw counts (features x samples) for species A and
#'   B. Sample names must be unique across both.
#' @param orthologs One-to-one ortholog table, e.g. from
#'   [reciprocal_best_hits()].
#' @param gene_a_col,gene_b_col Columns of `orthologs` with the gene ids.
#' @param id_col Column used as row name; if `NULL`, rows are named
#'   `"<gene_a>|<gene_b>"`.
#' @param lengths_a,lengths_b Named numeric vectors of lengths (bp). Required
#'   when `length_normalize = TRUE`.
#' @param length_normalize Apply `Q = round(E / L_kb)` per species, as in the
#'   original publication. This changes the count distribution assumed by
#'   DESeq2; keep `FALSE` unless replicating that analysis.
#' @return A count matrix (orthologs x samples of A followed by B).
#' @export
merge_ortholog_counts <- function(counts_a, counts_b, orthologs, gene_a_col = "gene_a", gene_b_col = "gene_b",
                                  id_col = "ortholog_id", lengths_a = NULL, lengths_b = NULL,
                                  length_normalize = FALSE) {
  shared <- intersect(colnames(counts_a), colnames(counts_b))
  if (length(shared)) stop("Sample names shared by both matrices: ", paste(utils::head(shared, 5), collapse = ", "))
  ga <- as.character(orthologs[[gene_a_col]])
  gb <- as.character(orthologs[[gene_b_col]])
  if (anyDuplicated(ga) || anyDuplicated(gb)) stop("The ortholog table must be one-to-one")
  present <- ga %in% rownames(counts_a) & gb %in% rownames(counts_b)
  if (any(!present)) warning(sum(!present), " ortholog pair(s) absent from the count matrices were dropped")
  rid <- if (is.null(id_col)) paste(ga, gb, sep = "|") else as.character(orthologs[[id_col]])
  ga <- ga[present]; gb <- gb[present]; rid <- rid[present]
  a <- as.matrix(counts_a[ga, , drop = FALSE]); storage.mode(a) <- "double"
  b <- as.matrix(counts_b[gb, , drop = FALSE]); storage.mode(b) <- "double"
  if (length_normalize) {
    if (is.null(lengths_a) || is.null(lengths_b)) stop("length_normalize = TRUE requires lengths_a and lengths_b")
    la <- lengths_a[ga] / 1000
    lb <- lengths_b[gb] / 1000
    if (anyNA(la) || anyNA(lb)) stop("Missing lengths for some orthologs")
    a <- round(a / la)
    b <- round(b / lb)
  }
  m <- cbind(a, b)
  rownames(m) <- rid
  if (all(abs(m - round(m)) < 1e-8)) storage.mode(m) <- "integer"
  m
}
