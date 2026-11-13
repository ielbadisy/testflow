test_that("testgroups recommends ANOVA and returns posthoc and effect size", {
  set.seed(4)
  dat <- .bt(y = rnorm(90), g = rep(c("a", "b", "c"), each = 30))
  x <- testgroups(dat, y, g)
  expect_s3_class(x, "testflow_groups")
  expect_equal(x$recommended$test, "One-way ANOVA")
  expect_true(!is.null(x$posthoc))
  expect_equal(x$posthoc$method, "Tukey HSD")
  expect_true(nrow(x$effect_size) >= 1)
  z <- testgroups(y ~ g, data = dat)
  expect_s3_class(z, "testflow_groups")
})

test_that("testgroups can recommend Kruskal", {
  dat <- .bt(y = c(rexp(25), rexp(25), rexp(25)), g = rep(c("a", "b", "c"), each = 25))
  x <- testgroups(dat, y, g)
  expect_equal(x$recommended$test, "Kruskal-Wallis test")
  expect_equal(x$posthoc$method, "Pairwise Wilcoxon rank-sum tests (BH-adjusted)")
})

test_that("testgroups selects Welch posthoc when variances differ", {
  set.seed(5)
  dat <- .bt(
    y = c(rnorm(50, 0, 1), rnorm(50, 0.5, 4), rnorm(50, 1, 7)),
    g = rep(c("a", "b", "c"), each = 50)
  )
  x <- testgroups(dat, y, g)
  expect_equal(x$recommended$test, "Welch ANOVA")
  expect_equal(x$posthoc$method, "Pairwise Welch t-tests (BH-adjusted)")
})
