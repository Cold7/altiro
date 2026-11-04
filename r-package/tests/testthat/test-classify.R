fixture_tables <- function() {
  keys <- names(contrast_keys())
  paths <- setNames(lapply(keys, function(k) {
    system.file("extdata", paste0("contrast_", k, ".tsv"), package = "orthoresponse")
  }), keys)
  read_contrasts(paths)
}

expected_labels <- function() {
  p <- system.file("extdata", "expected_labels.tsv", package = "orthoresponse")
  e <- utils::read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(e) <- e[[1]]
  e
}

test_that("labels match the shared fixture (same as the Python implementation)", {
  tabs <- fixture_tables()
  e <- expected_labels()
  expect_identical(classify_responses(tabs)[rownames(e), "label"], e$paper)
  expect_identical(classify_responses(tabs, preset = "strict")[rownames(e), "label"], e$strict)
  expect_identical(classify_responses(tabs, equivalence_margin = 0.5)[rownames(e), "label"],
                   e[["paper_equivalence_0.5"]])
})

one <- function(a, b, basal, mb, ma, ...) {
  vals <- list(a, b, basal, mb, ma)
  tabs <- setNames(lapply(vals, function(v) data.frame(padj = v[1], log2FoldChange = v[2], row.names = "g")),
                   names(contrast_keys()))
  classify_responses(tabs, ...)$label
}

test_that("single-feature rules follow the published decision tree", {
  ns <- 0.5; de <- 0.001
  expect_equal(one(c(ns, .1), c(de, 2), c(de, -2), c(ns, .1), c(de, 2)), "stress_ready_1")
  expect_equal(one(c(ns, .1), c(de, -2), c(de, 2), c(ns, .1), c(de, -2)), "stress_ready_2")
  expect_equal(one(c(de, 2), c(ns, .1), c(de, 2), c(de, -2), c(ns, .1)), "stress_ready_3")
  expect_equal(one(c(de, -2), c(ns, .1), c(de, -2), c(de, 2), c(ns, .1)), "stress_ready_4")
  expect_equal(one(c(de, 2), c(ns, .1), c(ns, .1), c(de, 2), c(de, 2)), "unique_1")
  expect_equal(one(c(ns, .1), c(de, -2), c(ns, .1), c(de, -2), c(ns, .1)), "unique_4")
  expect_equal(one(c(de, 2), c(de, 2), c(ns, 0), c(de, 2), c(de, 2)), "shared_1")
  expect_equal(one(c(de, -2), c(de, 2), c(ns, 0), c(ns, 0), c(ns, 0)), "opposite_2")
  expect_equal(one(c(ns, 0), c(ns, 0), c(de, 1), c(de, 1), c(de, -1)), "no_response_1")
  expect_equal(one(c(ns, 0), c(ns, 0), c(NA, 1), c(ns, 0), c(ns, 0)), "no_response_3")
})

test_that("strict preset requires a basal difference", {
  args <- list(c(0.001, 1), c(0.6, 0.1), c(0.3, 0.4), c(0.2, 0.5), c(0.4, 0.3))
  expect_equal(do.call(one, args), "stress_ready_3")
  expect_equal(do.call(one, c(args, preset = "strict")), "unique_1")
})

test_that("missing values are flagged or unclassified", {
  tabs <- fixture_tables()
  paper <- classify_responses(tabs)
  strict <- classify_responses(tabs, na_action = "unclassified")
  flagged <- grepl("na_treated_as_not_de", paper$flags)
  expect_true(any(flagged))
  expect_true(all(strict$model[flagged] == "unclassified"))
  expect_identical(strict$label[!flagged], paper$label[!flagged])
})

test_that("BH adjustment and TOST agree with closed forms", {
  res <- equivalence_test(c(0, 0.3, 1.2), c(0.1, 0.2, 0.2), 0.5)
  lfc <- c(0, 0.3, 1.2); se <- c(0.1, 0.2, 0.2)
  expect_equal(res$pvalue, pmax(pnorm((lfc + 0.5) / se, lower.tail = FALSE), pnorm((lfc - 0.5) / se)))
  expect_equal(stats::p.adjust(c(0.01, 0.04, NA, 0.03, 0.2, 0.5), "BH"),
               c(0.05, 0.0666666667, NA, 0.0666666667, 0.25, 0.5), tolerance = 1e-8)
})

test_that("summaries use explicit denominators", {
  cls <- classify_responses(fixture_tables(), group_names = c("sens", "tol"))
  s <- summarize_classification(cls)
  expect_equal(sum(s$n), nrow(cls))
  expect_equal(sum(s$pct_of_responsive, na.rm = TRUE), 100)
  a <- primed_asymmetry(cls)
  expect_equal(a$responding_group, c("sens", "tol"))
  expect_true(a$fisher_pvalue[1] >= 0 && a$fisher_pvalue[1] <= 1)
})

test_that("wide input gives the same result", {
  tabs <- fixture_tables()
  wide <- do.call(cbind, lapply(names(tabs), function(k) {
    d <- tabs[[k]][, c("log2FoldChange", "padj")]
    colnames(d) <- paste0(k, "_", colnames(d))
    d
  }))
  expect_identical(classify_responses(wide)$label, classify_responses(tabs)$label)
})

test_that("options are validated", {
  expect_error(classification_options("paper", na_action = "drop"))
  expect_error(classification_options("paper", alpha = 2))
})
