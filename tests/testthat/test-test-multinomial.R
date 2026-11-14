# Tests for testmultinomial() result layout.

test_that("testmultinomial results live in primary_test and per-category posthoc", {
  dat <- data.frame(g = factor(rep(c("a", "b", "c"), times = c(30, 20, 10))))
  x <- testmultinomial(dat, g, p = c(0.5, 0.3, 0.2), plot = FALSE)
  ref <- stats::chisq.test(c(30, 20, 10), p = c(0.5, 0.3, 0.2))
  expect_equal(x$primary_test$p.value, unname(ref$p.value))
  expect_equal(x$posthoc$level, c("a", "b", "c"))
  expect_equal(x$posthoc$p[1], stats::binom.test(30, 60, 0.5)$p.value)
})
