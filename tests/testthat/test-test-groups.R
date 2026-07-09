test_that("test_groups recommends ANOVA and returns posthoc and effect size", {
  set.seed(4)
  dat <- tibble::tibble(y = rnorm(90), g = rep(c("a", "b", "c"), each = 30))
  x <- test_groups(dat, y, g)
  expect_s3_class(x, "testflow_groups")
  expect_equal(x$recommended$test, "One-way ANOVA")
  expect_true(!is.null(x$posthoc))
  expect_equal(x$posthoc$method, "Tukey HSD")
  expect_true(nrow(x$effect_size) >= 1)
  z <- test_groups(y ~ g, data = dat)
  expect_s3_class(z, "testflow_groups")
})

test_that("test_groups can recommend Kruskal", {
  dat <- tibble::tibble(y = c(rexp(25), rexp(25), rexp(25)), g = rep(c("a", "b", "c"), each = 25))
  x <- test_groups(dat, y, g)
  expect_equal(x$recommended$test, "Kruskal-Wallis test")
  expect_equal(x$posthoc$method, "Pairwise Wilcoxon rank-sum tests (BH-adjusted)")
})

test_that("test_groups selects Welch posthoc when variances differ", {
  set.seed(5)
  dat <- tibble::tibble(
    y = c(rnorm(50, 0, 1), rnorm(50, 0.5, 4), rnorm(50, 1, 7)),
    g = rep(c("a", "b", "c"), each = 50)
  )
  x <- test_groups(dat, y, g)
  expect_equal(x$recommended$test, "Welch ANOVA")
  expect_equal(x$posthoc$method, "Pairwise Welch t-tests (BH-adjusted)")
})
