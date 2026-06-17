test_that("two groups recommend Welch or Wilcoxon appropriately", {
  set.seed(2)
  welch_dat <- tibble::tibble(y = c(rnorm(60, 0, 1), rnorm(60, 0, 4)), g = rep(c("a", "b"), each = 60))
  wilcox_dat <- tibble::tibble(y = c(rexp(35), rexp(35, rate = 0.5)), g = rep(c("a", "b"), each = 35))
  x <- test_two_groups(welch_dat, y, g, plot = FALSE)
  y <- test_two_groups(wilcox_dat, y, g, plot = FALSE)
  expect_equal(x$recommended$test, "Welch t-test")
  expect_equal(y$recommended$test, "Wilcoxon rank-sum test")
})
