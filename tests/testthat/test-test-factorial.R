test_that("test_factorial supports formula-first calls", {
  dat <- make_cardio_data(90)
  x <- test_factorial(sbp_3m ~ sex * treatment, data = dat)
  expect_s3_class(x, "testflow_factorial")
  expect_equal(x$outcome, "sbp_3m")
  expect_equal(x$recommended$test, "Factorial ANOVA")
  expect_true("null_hypothesis" %in% names(x$primary_test))
  expect_true("effect_size" %in% names(x$primary_test))
  expect_s3_class(plot(x), "ggplot")
})

test_that("test_factorial's type parameter actually selects Type I/II/III sums of squares", {
  skip_if_not_installed("car")
  set.seed(1)
  n <- c(10, 15, 8, 20)
  dat <- data.frame(
    y = c(rnorm(10, 10, 2), rnorm(15, 12, 2), rnorm(8, 9, 2), rnorm(20, 14, 2)),
    A = factor(rep(c("a1", "a1", "a2", "a2"), n)),
    B = factor(rep(c("b1", "b2", "b1", "b2"), n))
  )
  lm_fit2 <- lm(y ~ A * B, data = dat)
  ref2 <- car::Anova(lm_fit2, type = 2)
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  lm_fit3 <- lm(y ~ A * B, data = dat)
  ref3 <- car::Anova(lm_fit3, type = 3)
  options(old)

  x1 <- test_factorial(y ~ A * B, data = dat, type = 1, plot = FALSE)
  x2 <- test_factorial(y ~ A * B, data = dat, type = 2, plot = FALSE)
  x3 <- test_factorial(y ~ A * B, data = dat, type = 3, plot = FALSE)

  a_p <- function(x) x$alternative_tests$anova_table$p.value[x$alternative_tests$anova_table$term == "A"]
  expect_equal(a_p(x2), ref2["A", "Pr(>F)"])
  expect_equal(a_p(x3), ref3["A", "Pr(>F)"])
  expect_false(isTRUE(all.equal(a_p(x1), a_p(x2))))
  expect_false(isTRUE(all.equal(a_p(x2), a_p(x3))))

  # eta squared must be computed from the same SS type as the reported test,
  # not always Type I regardless of `type`.
  eta_a <- function(x) x$effect_size$estimate[x$effect_size$term == "A"]
  expect_false(isTRUE(all.equal(eta_a(x1), eta_a(x2))))
})
