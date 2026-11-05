test_that("plot_model_patterns accepts custom condition labels and title", {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(plot_model_patterns())
  expect_null(plot_model_patterns(condition_labels = c("Ctrl", "Drought"), group_names = c("Sly", "Spenn"),
                                  title = "Tomato drought"))
  expect_error(plot_model_patterns(condition_labels = "Control"), "exactly two")
  expect_null(plot_model_patterns(model_labels = c(primed = "Pre-adapted", no_response = "Unchanged")))
  expect_null(plot_model_patterns(model_labels = c(stress_ready = "Stress-ready")))   # legacy key
  expect_null(plot_model_patterns(model_labels = list(opposite = "Inverse")))
  expect_error(plot_model_patterns(model_labels = c(stressready = "x")), "model_labels")
  expect_error(plot_model_patterns(model_labels = "Primed"), "model_labels")
})

test_that("save_plots writes the figures", {
  keys <- names(contrast_keys())
  paths <- setNames(lapply(keys, function(k) system.file("extdata", paste0("contrast_", k, ".tsv"),
                                                           package = "orthoresponse")), keys)
  cls <- classify_responses(read_contrasts(paths))
  out <- save_plots(cls, tempfile(), formats = "png", condition_labels = c("C", "D"))
  expect_equal(basename(out), c("model_patterns.png", "model_summary.png"))
  expect_true(all(file.size(out) > 1000))
})
