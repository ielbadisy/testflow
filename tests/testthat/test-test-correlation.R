test_that("test_correlation returns a testflow object", {
  dat <- make_cardio_data(60)
  x <- test_correlation(sbp_3m ~ age, data = dat)
  expect_s3_class(x, "testflow_correlation")
  expect_s3_class(plot(x), "ggplot")
})

test_that("test_correlation_matrix uses genuinely pairwise-complete observations", {
  set.seed(1)
  d <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  d$a[1:5] <- NA
  d$c[26:30] <- NA

  x <- suppressWarnings(test_correlation_matrix(d, c(a, b, c), method = "pearson", plot = FALSE))
  bc <- x$alternative_tests$correlation_table
  bc_row <- bc[bc$var1 == "b" & bc$var2 == "c", ]

  expect_equal(bc_row$estimate, cor(d$b, d$c, use = "pairwise.complete.obs"))
  expect_false(isTRUE(all.equal(
    bc_row$estimate,
    cor(na.omit(d)$b, na.omit(d)$c)
  )))
})
