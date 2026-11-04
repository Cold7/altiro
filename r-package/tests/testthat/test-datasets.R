test_that("the rice example loads with matching design", {
  ex <- load_example("rice_salt")
  expect_equal(dim(ex$counts), c(36844L, 12L))
  expect_identical(colnames(ex$counts), rownames(ex$samples))
  expect_equal(sum(as.numeric(ex$counts)), 198199444)
  d <- build_design(ex$samples, ex$params$group_col, ex$params$condition_col, ex$params$group_a,
                    ex$params$group_b, ex$params$control, ex$params$stress)
  expect_equal(as.vector(table(d$cell)), c(3L, 3L, 3L, 3L))
  expect_true("rice_salt" %in% list_examples()$name)
  expect_error(load_example("nope"), "Unknown example")
})

test_that("export_example writes the files", {
  out <- export_example(tempfile())
  expect_true(all(file.exists(out)))
})

test_that("the ovarian example loads with matching design", {
  ex <- load_example("ovarian_brequinar")
  expect_equal(dim(ex$counts), c(59453L, 12L))
  expect_identical(colnames(ex$counts), rownames(ex$samples))
  expect_equal(sum(as.numeric(ex$counts)), 281005562)
  d <- build_design(ex$samples, ex$params$group_col, ex$params$condition_col, ex$params$group_a,
                    ex$params$group_b, ex$params$control, ex$params$stress)
  expect_equal(as.vector(table(d$cell)), c(3L, 3L, 3L, 3L))
  expect_setequal(list_examples()$name, c("rice_salt", "ovarian_brequinar"))
})
