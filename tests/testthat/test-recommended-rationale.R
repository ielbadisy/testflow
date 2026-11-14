# Tests that auto-chosen tests carry a specific rationale.

test_that("auto-choice workflows explain the recommended test with assumption p-values", {
  set.seed(3)
  dat <- data.frame(y = c(rnorm(15), rexp(15) * 5, rnorm(15, 1)), g = factor(rep(c("a", "b", "c"), each = 15)))
  x <- testgroups(dat, y, g, plot = FALSE)
  expect_type(x$recommended$rationale, "character")
  expect_match(x$recommended$rationale, "Shapiro-Wilk", fixed = TRUE)
  expect_match(x$recommended$rationale, "Levene", fixed = TRUE)
  expect_match(x$recommended$rationale, x$recommended$test, fixed = TRUE)

  d2 <- dat[dat$g != "c", ]
  d2$g <- droplevels(d2$g)
  y <- testtwogroups(d2, y, g, plot = FALSE)
  expect_match(y$recommended$rationale, "Shapiro-Wilk", fixed = TRUE)
  expect_match(y$recommended$rationale, "Levene", fixed = TRUE)

  z <- testonesample(data.frame(v = rnorm(20)), v, mu = 0, plot = FALSE)
  expect_match(z$recommended$rationale, "Shapiro-Wilk", fixed = TRUE)

  p <- testpaired(data.frame(b = rnorm(20), a = rnorm(20)), b, a, plot = FALSE)
  expect_match(p$recommended$rationale, "differences", fixed = TRUE)

  cc <- testcorrelation(data.frame(u = rnorm(30), w = rnorm(30)), u, w, plot = FALSE)
  expect_match(cc$recommended$rationale, "Shapiro-Wilk", fixed = TRUE)

  long <- data.frame(id = factor(rep(1:12, 3)), time = factor(rep(c("t1", "t2", "t3"), each = 12)), y = rnorm(36))
  r <- testrepeatedlong(long, y, time, id, plot = FALSE)
  expect_match(r$recommended$rationale, "Shapiro-Wilk", fixed = TRUE)
})

test_that("print() shows the recommendation rationale", {
  x <- testonesample(data.frame(v = c(1.2, 0.4, -0.3, 2.1, 0.8, -1.1, 0.5, 1.7, 0.2, -0.6)), v, plot = FALSE)
  out <- paste(capture.output(print(x)), collapse = "\n")
  expect_match(out, x$recommended$rationale, fixed = TRUE)
})
