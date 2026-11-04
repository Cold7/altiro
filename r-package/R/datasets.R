.examples <- list(
  rice_salt = list(
    counts = "rice_GSE167342_counts.tsv.gz",
    samples = "rice_GSE167342_samples.tsv",
    params = list(group_col = "line", condition_col = "treatment", group_a = "PL12", group_b = "IL22",
                  control = "Control", stress = "Salt"),
    description = paste("Rice leaves (GEO GSE167342): salt-sensitive PL12 (A) vs salt-tolerant Saltol",
                        "introgression line IL22 (B), control vs 24 h 80 mM NaCl, 3 replicates.",
                        "Raw counts, 36,844 genes x 12 samples.")
  ),
  ovarian_brequinar = list(
    counts = "ovarian_GSE283381_counts.tsv.gz",
    samples = "ovarian_GSE283381_samples.tsv",
    params = list(group_col = "line", condition_col = "treatment", group_a = "CisS", group_b = "CisR",
                  control = "Vehicle", stress = "Brequinar"),
    description = paste("Human ovarian cancer cells (GEO GSE283381): OVCAR5 cisplatin-sensitive (A) vs",
                        "cisplatin-resistant (B), vehicle vs the DHODH inhibitor brequinar, 3 replicates.",
                        "Raw HTSeq counts, 59,453 genes x 12 samples.")
  )
)

.get_example <- function(name) {
  if (!name %in% names(.examples)) {
    stop("Unknown example '", name, "'; available: ", paste(names(.examples), collapse = ", "))
  }
  .examples[[name]]
}

#' Example data sets shipped with the package
#'
#' `rice_salt`: rice leaves, GEO GSE167342. Salt-sensitive *japonica* cultivar
#' PL12 (group A) and the salt-tolerant introgression line IL22 (group B,
#' derived from FL478 x PL12, carries the *Saltol* QTL), control vs 24 h
#' 80 mM NaCl, 3 replicates per group. Raw gene-level counts (36,844 genes x
#' 12 samples) as deposited by the authors in GEO
#' (`GSE167342_Raw_counts.xlsx`); columns S40-S72 are the authors' library
#' names, described in the sample table.
#'
#' `ovarian_brequinar`: human ovarian cancer cells, GEO GSE283381. OVCAR5
#' cisplatin-sensitive cells (CisS, group A) and their cisplatin-resistant
#' derivative (CisR, group B), treated with vehicle or the DHODH inhibitor
#' brequinar, 3 replicates per group. Raw HTSeq counts deposited by the
#' authors (one file per sample in `GSE283381_RAW.tar`), merged into one
#' matrix (59,453 genes x 12 samples); columns are the authors' library names.
#'
#' @return `list_examples()`: a data frame with names and descriptions.
#' @examples
#' list_examples()
#' @export
list_examples <- function() {
  data.frame(name = names(.examples), description = vapply(.examples, `[[`, "", "description"),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' @rdname list_examples
#' @param name Name of the example (see `list_examples()`).
#' @return `example_files()`: named character vector with the paths of the
#'   `counts` and `samples` files.
#' @export
example_files <- function(name = "rice_salt") {
  info <- .get_example(name)
  c(counts = system.file("extdata", info$counts, package = "orthoresponse", mustWork = TRUE),
    samples = system.file("extdata", info$samples, package = "orthoresponse", mustWork = TRUE))
}

#' @rdname list_examples
#' @return `load_example()`: a list with `counts` (integer matrix, genes x
#'   samples), `samples` (data frame), `params` (the matching arguments for
#'   [orthoresponse()]) and `description`.
#' @examples
#' ex <- load_example("rice_salt")
#' dim(ex$counts)
#' ex$samples
#' \dontrun{
#' res <- do.call(orthoresponse, c(list(ex$counts, ex$samples), ex$params))
#' res
#' }
#' @export
load_example <- function(name = "rice_salt") {
  info <- .get_example(name)
  f <- example_files(name)
  counts <- utils::read.delim(f[["counts"]], row.names = 1, check.names = FALSE)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"
  samples <- utils::read.delim(f[["samples"]], row.names = 1, stringsAsFactors = FALSE)
  list(name = name, counts = counts, samples = samples, params = info$params, description = info$description)
}

#' @rdname list_examples
#' @param outdir Folder where the files are copied (`counts.tsv.gz`,
#'   `samples.tsv`); useful as a template for your own input.
#' @return `export_example()`: invisibly, the paths of the copied files.
#' @export
export_example <- function(outdir = ".", name = "rice_salt") {
  f <- example_files(name)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  out <- c(counts = file.path(outdir, "counts.tsv.gz"), samples = file.path(outdir, "samples.tsv"))
  ok <- file.copy(f, out, overwrite = TRUE)
  if (!all(ok)) stop("Could not copy the example files to ", outdir)
  invisible(out)
}
