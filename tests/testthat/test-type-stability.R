# Regression tests for testflow 1.0.4: integer outcomes with odd and even
# group sizes, and large samples (n > 5000).

test_that("test_groups() handles integer outcomes with odd and even group sizes", {
  set.seed(1)
  d <- data.frame(
    y = c(rpois(101, 3), rpois(100, 3), rpois(51, 3)),
    g = factor(rep(c("a", "b", "c"), c(101, 100, 51)))
  )
  expect_true(is.integer(d$y))
  res <- suppressWarnings(test_groups(y ~ g, data = d, plot = FALSE))
  expect_s3_class(res, "testflow")
  expect_type(res$descriptives$median, "double")
})

test_that("sumtab() handles an integer variable with an all-missing group", {
  d <- data.frame(y = c(1:5, rep(NA_integer_, 4)), g = factor(rep(c("a", "b"), c(5, 4))))
  expect_no_error(suppressWarnings(sumtab(~ y | g, data = d)))
})

test_that(".rbind_rows() widens integer/double and all-NA columns only", {
  out <- testflow:::.rbind_rows(data.frame(a = 1L, b = NA), data.frame(a = 1.5, b = "x"))
  expect_type(out$a, "double")
  expect_type(out$b, "character")
  expect_error(testflow:::.rbind_rows(data.frame(a = 1L), data.frame(a = "z")))
})

test_that("regression and factorial workflows run with n > 5000", {
  set.seed(2)
  n <- 6000
  d <- data.frame(
    y = rnorm(n), x = rnorm(n),
    g = factor(sample(c("a", "b", "c"), n, TRUE)),
    h = factor(sample(c("u", "v"), n, TRUE))
  )
  lin <- suppressWarnings(test_linear_regression(y ~ x + g, data = d, plot = FALSE))
  fac <- suppressWarnings(test_factorial(y ~ g * h, data = d, plot = FALSE))
  expect_s3_class(lin, "testflow")
  expect_s3_class(fac, "testflow")
  expect_true(any(lin$assumptions$status == "not checked"))
})

test_that("test_categorical() falls back to a Monte Carlo Fisher p-value on large tables", {
  set.seed(3)
  n <- 12500
  d <- data.frame(
    x = factor(sample(c("no", "yes"), n, TRUE)),
    g = factor(sample(c("a", "b", "c"), n, TRUE))
  )
  res <- suppressWarnings(test_categorical(x ~ g, data = d, plot = FALSE))
  expect_s3_class(res, "testflow")
  expect_true(is.finite(res$alternative_tests$fisher$p.value))
})
