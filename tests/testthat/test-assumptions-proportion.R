test_that("proportion workflow chooses exact binomial when expected counts are small", {
  dat <- tibble::tibble(x = c("yes", "yes", "no", "no", "no"))
  x <- test_proportion(dat, x, success = "yes", p = 0.8, plot = FALSE)
  expect_equal(x$recommended$test, "Exact binomial test")
})
