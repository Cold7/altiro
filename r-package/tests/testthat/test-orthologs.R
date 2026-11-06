blast <- function(rows) {
  d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(d) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  for (k in 3:12) d[[k]] <- as.numeric(d[[k]])
  d
}

test_that("reciprocal best hits are found", {
  ab <- blast(list(
    c("a1.1", "b1.1", 90, 100, 1, 0, 1, 100, 1, 100, 1e-50, 300),
    c("a1.1", "b2.1", 60, 100, 1, 0, 1, 100, 1, 100, 1e-20, 120),
    c("a2.1", "b2.1", 80, 100, 1, 0, 1, 100, 1, 100, 1e-40, 250),
    c("a3.1", "b2.1", 70, 100, 1, 0, 1, 100, 1, 100, 1e-30, 200)
  ))
  ba <- blast(list(
    c("b1.1", "a1.1", 90, 100, 1, 0, 1, 100, 1, 100, 1e-50, 300),
    c("b2.1", "a2.1", 80, 100, 1, 0, 1, 100, 1, 100, 1e-40, 250),
    c("b2.1", "a3.1", 70, 100, 1, 0, 1, 100, 1, 100, 1e-30, 200)
  ))
  f_ab <- tempfile(); f_ba <- tempfile()
  utils::write.table(ab, f_ab, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  utils::write.table(ba, f_ba, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  r <- reciprocal_best_hits(f_ab, f_ba, strip_pattern = "\\.\\d+$")
  expect_equal(r$gene_a, c("a1", "a2"))
  expect_equal(r$gene_b, c("b1", "b2"))
  expect_equal(r$ortholog_id, c("ortholog_1", "ortholog_2"))
})

test_that("best_hits breaks ties by subject id", {
  d <- blast(list(c("q", "s2", 1, 1, 1, 1, 1, 1, 1, 1, 1e-5, 50), c("q", "s1", 1, 1, 1, 1, 1, 1, 1, 1, 1e-5, 50)))
  expect_equal(best_hits(d)$sseqid, "s1")
})

test_that("merge and length normalisation", {
  ca <- data.frame(a_c1 = c(100, 10), a_s1 = c(200, 20), row.names = c("a1", "a2"))
  cb <- data.frame(b_c1 = c(50, 5), b_s1 = c(60, 6), row.names = c("b1", "b2"))
  orth <- data.frame(ortholog_id = c("o1", "o2"), gene_a = c("a1", "a2"), gene_b = c("b1", "b2"))
  m <- merge_ortholog_counts(ca, cb, orth)
  expect_equal(dim(m), c(2L, 4L))
  expect_equal(m["o1", "b_s1"], 60L)
  m2 <- merge_ortholog_counts(ca, cb, orth, lengths_a = c(a1 = 2000, a2 = 500), lengths_b = c(b1 = 1000, b2 = 4000),
                              length_normalize = TRUE)
  expect_equal(unname(m2[, "a_c1"]), c(50L, 20L))
  expect_equal(m2["o2", "b_c1"], 1L)
  # normalized values stay integer counts (rounded), as DESeq2 requires
  expect_true(is.integer(m2))
  m3 <- merge_ortholog_counts(ca, cb, orth, lengths_a = c(a1 = 3000, a2 = 700),
                              lengths_b = c(b1 = 1300, b2 = 4100), length_normalize = TRUE)
  expect_true(is.integer(m3))
  orth2 <- data.frame(ortholog_id = c("o1", "o2"), gene_a = c("a1", "a2"), gene_b = c("b1", "b1"))
  expect_error(merge_ortholog_counts(ca, cb, orth2))
})
