test_that("Spearman correlation does not require normality", {
  dat <- make_assumption_data(60)
  x <- test_correlation(sbp_3m ~ age, data = dat, method = "spearman", plot = FALSE)
  expect_true(any(grepl("Monotonic relationship", x$assumptions$name)))
  expect_true(any(grepl("Normality", x$assumptions$name)))
  expect_true(any(x$assumptions$status == "not required"))
})

test_that("Pearson correlation reports linearity normality and outliers", {
  dat <- make_assumption_data(60)
  x <- test_correlation(sbp_3m ~ age, data = dat, method = "pearson", plot = FALSE)
  expect_true(any(grepl("Linearity", x$assumptions$name)))
  expect_true(any(x$assumptions$variable == "age"))
  expect_true(any(x$assumptions$variable == "sbp_3m"))
})
