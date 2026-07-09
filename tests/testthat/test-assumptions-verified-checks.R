test_that("check_symmetry matches the closed-form Cabilio-Masaro statistic", {
  set.seed(11)
  x <- rexp(30, 1)
  n <- length(x)
  expected_stat <- sqrt(n) * (mean(x) - median(x)) / (sd(x) * sqrt(pi / 2 - 1))
  expected_p <- 2 * (1 - pnorm(abs(expected_stat)))

  out <- check_symmetry(x, alpha = 0.05)
  expect_equal(out$statistic, expected_stat)
  expect_equal(out$p_value, expected_p)
  expect_identical(out$status, ifelse(expected_p >= 0.05, "acceptable", "not acceptable"))
})

test_that("check_symmetry flags a symmetric sample as acceptable", {
  set.seed(12)
  x <- rnorm(200, 0, 1)
  out <- check_symmetry(x, alpha = 0.05)
  expect_equal(out$status, "acceptable")
})

test_that("check_linearity matches the quadratic-term F-test from anova()", {
  set.seed(13)
  x <- runif(40, 0, 10)
  y <- 2 + 0.5 * x + 0.3 * x^2 + rnorm(40, 0, 3)
  fit_linear <- lm(y ~ x)
  fit_quad <- lm(y ~ x + I(x^2))
  cmp <- anova(fit_linear, fit_quad)

  out <- check_linearity(x, y, alpha = 0.05)
  expect_equal(out$statistic, cmp$F[2])
  expect_equal(out$p_value, cmp[["Pr(>F)"]][2])
  expect_equal(out$status, "warning")
})

test_that("check_linearity accepts a genuinely linear relationship", {
  set.seed(14)
  x <- runif(60, 0, 10)
  y <- 3 + 2 * x + rnorm(60, 0, 1)
  out <- check_linearity(x, y, alpha = 0.05)
  expect_equal(out$status, "acceptable")
})

test_that("check_sphericity matches base R mauchly.test() for 3+ conditions", {
  set.seed(15)
  n <- 20
  t1 <- rnorm(n, 10, 2)
  t2 <- t1 + rnorm(n, 1, 1.5)
  t3 <- t2 + rnorm(n, 1, 3)
  wide <- cbind(t1, t2, t3)
  mlm_fit <- lm(wide ~ 1)
  mt <- mauchly.test(mlm_fit, X = ~1)

  out <- check_sphericity(wide, alpha = 0.05)
  expect_equal(out$statistic, unname(mt$statistic))
  expect_equal(out$p_value, unname(mt$p.value))
})

test_that("check_sphericity reports not applicable for exactly two conditions", {
  wide <- cbind(a = rnorm(10), b = rnorm(10))
  out <- check_sphericity(wide, alpha = 0.05)
  expect_equal(out$status, "not applicable")
})

test_that("test_correlation_matrix applies BH correction and reports it", {
  set.seed(16)
  df <- data.frame(a = rnorm(30), b = rnorm(30), c = rnorm(30))
  x <- test_correlation_matrix(df, c(a, b, c), method = "pearson", plot = FALSE)
  tab <- x$alternative_tests$correlation_table
  expect_true("p.adj" %in% names(tab))
  expect_equal(tab$p.adj, p.adjust(tab$p, method = "BH"))
  expect_true(any(grepl("Multiple testing correction", x$assumptions$name)))
  expect_equal(x$assumptions$status[x$assumptions$name == "Multiple testing correction"], "acceptable")
})
