# Tests for testcorrelationmatrix() pairwise results.

test_that("testcorrelationmatrix exposes the pairwise table as posthoc", {
  dat <- mtcars[, c("mpg", "hp", "wt")]
  dat$hp[3] <- NA
  x <- suppressWarnings(testcorrelationmatrix(dat, c(mpg, hp, wt), plot = FALSE))
  expect_false(is.null(x$posthoc))
  expect_equal(nrow(x$posthoc), 3)
  expect_true(all(c("var1", "var2", "estimate", "p", "p.adj", "p.adjust.method", "n") %in% names(x$posthoc)))
  expect_equal(x$posthoc$estimate, x$alternative_tests$correlation_table$estimate)
  expect_equal(x$posthoc$n, c(31L, 32L, 31L))
  expect_equal(x$posthoc$p.adj, stats::p.adjust(x$posthoc$p, "BH"))
})
