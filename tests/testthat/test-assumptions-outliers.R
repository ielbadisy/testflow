test_that("outlier methods report different assumptions", {
  dat <- make_assumption_data(60)
  expect_warning(
    iqr <- test_outliers(c(sbp_3m), data = dat, method = "iqr", plot = FALSE),
    "screening workflow"
  )
  expect_warning(
    mahal <- test_outliers(c(sbp_3m, ldl, crp), data = dat, method = "mahalanobis", plot = FALSE),
    "screening workflow"
  )
  expect_true(any(grepl("Skewness sensitivity", iqr$assumptions$name)))
  expect_true(any(grepl("Approximate multivariate normality", mahal$assumptions$name)))
})
