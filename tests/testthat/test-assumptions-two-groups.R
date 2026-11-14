test_that("two groups recommend Welch or Wilcoxon appropriately", {
  set.seed(2)
  welch_dat <- .bt(y = c(rnorm(60, 0, 1), rnorm(60, 0, 4)), g = rep(c("a", "b"), each = 60))
  wilcox_dat <- .bt(y = c(rexp(35), rexp(35, rate = 0.5)), g = rep(c("a", "b"), each = 35))
  x <- testtwogroups(welch_dat, y, g, plot = FALSE)
  y <- testtwogroups(wilcox_dat, y, g, plot = FALSE)
  expect_equal(x$recommended$test, "Welch t-test")
  expect_equal(y$recommended$test, "Wilcoxon rank-sum test")
})

test_that("two-group variance checks name their methods; Levene drives the choice", {
  set.seed(2)
  d <- data.frame(y = rnorm(40), g = rep(c("a", "b"), 20))
  r <- testtwogroups(y ~ g, d, plot = FALSE)
  a <- r$assumptions
  expect_equal(a$method[a$name == "Variance homogeneity"], "Levene test")
  expect_equal(a$method[a$name == "Variance ratio check"], "F-test (var.test)")
})
