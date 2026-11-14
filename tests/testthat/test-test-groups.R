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

test_that("Welch ANOVA and one-way ANOVA results show both degrees of freedom", {
  set.seed(1)
  dat <- data.frame(y = c(rnorm(20, 0, 1), rnorm(20, 1, 3), rnorm(20, 2, 6)), g = rep(c("a", "b", "c"), each = 20))
  w <- testgroups(y ~ g, data = dat, plot = FALSE)
  expect_equal(w$recommended$test, "Welch ANOVA")
  ref <- oneway.test(y ~ g, data = dat)
  expect_equal(w$primary_test$num.df, unname(ref$parameter[1]))
  expect_equal(w$primary_test$den.df, unname(ref$parameter[2]))
  df_txt <- sprintf("df = %.2f, %.2f", ref$parameter[1], ref$parameter[2])
  expect_match(format_primary_result(w), df_txt, fixed = TRUE)
  expect_match(report(w), df_txt, fixed = TRUE)
  set.seed(4)
  a <- testgroups(y ~ g, data = data.frame(y = rnorm(90), g = rep(c("a", "b", "c"), each = 30)), plot = FALSE)
  expect_equal(a$recommended$test, "One-way ANOVA")
  expect_match(format_primary_result(a), "df = 2.00, 87.00", fixed = TRUE)
})
