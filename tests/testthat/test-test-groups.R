test_that("test_groups recommends ANOVA and returns posthoc and effect size", {
  set.seed(4)
  dat <- tibble::tibble(y = rnorm(90), g = rep(c("a", "b", "c"), each = 30))
  x <- test_groups(dat, y, g)
  expect_s3_class(x, "testflow_groups")
  expect_equal(x$recommended$test, "One-way ANOVA")
  expect_true(!is.null(x$posthoc))
  expect_true(nrow(x$effect_size) >= 1)
})

test_that("test_groups can recommend Kruskal", {
  dat <- tibble::tibble(y = c(rexp(25), rexp(25), rexp(25)), g = rep(c("a", "b", "c"), each = 25))
  x <- test_groups(dat, y, g)
  expect_equal(x$recommended$test, "Kruskal-Wallis test")
})
