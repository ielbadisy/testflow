test_that("outlier methods report different assumptions", {
  dat <- make_assumption_data(60)
  iqr <- test_outliers(dat, c(sbp_3m), method = "iqr", plot = FALSE)
  mahal <- test_outliers(dat, c(sbp_3m, ldl, crp), method = "mahalanobis", plot = FALSE)
  expect_true(any(grepl("Skewness sensitivity", iqr$assumptions$name)))
  expect_true(any(grepl("Approximate multivariate normality", mahal$assumptions$name)))
})
