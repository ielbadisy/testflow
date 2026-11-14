test_that("testfactorial supports formula-first calls", {
  dat <- makecardiodata(90)
  x <- testfactorial(sbp_3m ~ sex * treatment, data = dat)
  expect_s3_class(x, "testflow_factorial")
  expect_equal(x$outcome, "sbp_3m")
  expect_equal(x$recommended$test, "Factorial ANOVA")
  expect_true("null_hypothesis" %in% names(x$primary_test))
  expect_true("effect_size" %in% names(x$primary_test))
  expect_s3_class(plot(x), "ggplot")
})

test_that("testfactorial's type parameter actually selects Type I/II/III sums of squares", {
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

  x1 <- testfactorial(y ~ A * B, data = dat, type = 1, plot = FALSE)
  x2 <- testfactorial(y ~ A * B, data = dat, type = 2, plot = FALSE)
  x3 <- testfactorial(y ~ A * B, data = dat, type = 3, plot = FALSE)

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

test_that("testfactorial reports every term and checks variance across all cells", {
  skip_if_not_installed("car")
  dat <- data.frame(mpg = mtcars$mpg, cyl = factor(mtcars$cyl), am = factor(mtcars$am))
  x <- testfactorial(mpg ~ cyl * am, data = dat, plot = FALSE)
  ref <- car::Anova(lm(mpg ~ cyl * am, data = dat), type = 2)
  terms <- c("cyl", "am", "cyl:am")
  expect_equal(x$primary_test$term, terms)
  expect_equal(x$primary_test$statistic, unname(ref[terms, "F value"]))
  expect_equal(x$primary_test$p.value, unname(ref[terms, "Pr(>F)"]))
  expect_equal(x$primary_test$parameter, unname(ref[terms, "Df"]))
  expect_equal(x$primary_test$den.df, rep(ref["Residuals", "Df"], 3))
  expect_equal(x$primary_test$effect_size, x$effect_size$estimate)
  expect_match(x$primary_test$null_hypothesis[3], "interaction", fixed = TRUE)
  out <- paste(capture.output(print(x)), collapse = "\n")
  for (tm in terms) expect_match(out, paste0(tm, ": "), fixed = TRUE)
  expect_match(report(x), "cyl:am", fixed = TRUE)
  lev <- car::leveneTest(mpg ~ interaction(cyl, am), data = dat)
  a <- x$assumptions[x$assumptions$name == "Variance homogeneity", ]
  expect_equal(a$statistic, lev[["F value"]][1])
})
