simulate_counts <- function(seed = 3) {
  set.seed(seed)
  n <- 1500
  base <- exp(stats::runif(n, log(30), log(3000)))
  f <- matrix(1, n, 4)
  f[1:150, 2:4] <- 4          # stress_ready_3
  f[151:300, 2] <- 4          # unique_1
  f[301:450, c(2, 4)] <- 0.25 # shared_2
  f[451:600, 3:4] <- 3        # no_response_1
  f[601:675, 4] <- 4          # unique_3
  f[676:750, c(2, 4)] <- 4    # shared_1
  mu <- (base * f)[, rep(1:4, each = 3)]
  counts <- matrix(stats::rnbinom(length(mu), mu = mu, size = 1 / 0.03), nrow = n)
  samples <- paste(rep(c("sens", "tol"), each = 6), rep(rep(c("ctrl", "salt"), each = 3), 2), 1:3, sep = "_")
  dimnames(counts) <- list(paste0("gene", seq_len(n)), samples)
  meta <- data.frame(line = sub("_.*", "", samples), treatment = sub("^[^_]+_([^_]+)_.*", "\\1", samples),
                     row.names = samples)
  list(counts = counts, meta = meta)
}

test_that("design cells are built correctly", {
  s <- simulate_counts()
  d <- build_design(s$meta, "line", "treatment", "sens", "tol", "ctrl", "salt")
  expect_equal(as.vector(table(d$cell)), c(3L, 3L, 3L, 3L))
})

test_that("the pipeline recovers simulated models", {
  skip_if_not_installed("DESeq2")
  s <- simulate_counts()
  res <- suppressWarnings(orthoresponse(s$counts, s$meta, group_col = "line", condition_col = "treatment",
                                        group_a = "sens", group_b = "tol", control = "ctrl", stress = "salt"))
  lab <- res$classification[paste0("gene", 1:600), "label"]
  expect_gt(mean(lab[1:150] == "stress_ready_3"), 0.9)
  expect_gt(mean(lab[151:300] == "unique_1"), 0.9)
  expect_gt(mean(lab[301:450] == "shared_2"), 0.9)
  expect_gt(mean(lab[451:600] == "no_response_1"), 0.85)
  out <- tempfile()
  save_orthoresponse(res, out)
  expect_true(file.exists(file.path(out, "classification.tsv")))
  expect_true(file.exists(file.path(out, "contrasts", "contrast_basal.tsv")))
})
