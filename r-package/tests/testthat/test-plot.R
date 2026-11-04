test_that("plot_model_patterns accepts custom condition labels and title", {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(plot_model_patterns())
  expect_null(plot_model_patterns(condition_labels = c("Ctrl", "Drought"), group_names = c("Sly", "Spenn"),
                                  title = "Tomato drought"))
  expect_error(plot_model_patterns(condition_labels = "Control"), "exactly two")
})
