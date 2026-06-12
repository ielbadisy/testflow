test_that("plot.testflow returns stored plot", {
  dat <- make_cardio_data(80)
  x <- test_groups(dat, sbp_3m, treatment)
  expect_s3_class(plot(x), "ggplot")
})
