test_that("test_correlation returns a testflow object", {
  dat <- make_cardio_data(60)
  x <- test_correlation(dat, age, sbp_3m)
  expect_s3_class(x, "testflow_correlation")
  expect_s3_class(plot(x), "ggplot")
})
