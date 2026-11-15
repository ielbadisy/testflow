test_that("TOST for means equals two one-sided t.test() calls (Welch and Student)", {
  set.seed(2)
  d <- data.frame(y = c(rnorm(40, 10, 3), rnorm(50, 10.5, 4)), g = rep(c("a", "b"), c(40, 50)))
  for (ve in c(FALSE, TRUE)) {
    r <- testequivalence(y ~ g, d, lower = -2, upper = 2, var.equal = ve, plot = FALSE)
    pl <- t.test(y ~ g, d, mu = -2, alternative = "greater", var.equal = ve)
    pu <- t.test(y ~ g, d, mu = 2, alternative = "less", var.equal = ve)
    expect_equal(r$alternative_tests$one_sided$p.value, c(pl$p.value, pu$p.value))
    expect_equal(r$primary_test$p.value, max(pl$p.value, pu$p.value))
    ci90 <- t.test(y ~ g, d, var.equal = ve, conf.level = 0.90)$conf.int
    expect_equal(c(r$primary_test$conf.low, r$primary_test$conf.high), as.numeric(ci90))
    expect_equal(unname(r$primary_test$parameter), unname(pl$parameter))
  }
  expect_output(print(r), "90% CI")
})

test_that("non-inferiority uses one bound only", {
  set.seed(4)
  d <- data.frame(y = c(rnorm(30, 5), rnorm(30, 5.2)), g = rep(c("new", "old"), each = 30))
  r <- testequivalence(y ~ g, d, lower = -0.8, upper = Inf, plot = FALSE)
  ref <- t.test(y ~ g, d, mu = -0.8, alternative = "greater")
  expect_equal(r$primary_test$p.value, ref$p.value)
  expect_equal(r$primary_test$conf.low, ref$conf.int[1])
  expect_equal(r$objective, "noninferiority")
})

test_that("TOST for proportions matches the Wald formulas by hand", {
  set.seed(5)
  d <- data.frame(cure = c(rbinom(120, 1, 0.80), rbinom(110, 1, 0.78)), arm = rep(c("A", "B"), c(120, 110)))
  r <- testequivalence(cure ~ arm, d, lower = -0.1, upper = 0.1, success = 1, plot = FALSE)
  p1 <- mean(d$cure[d$arm == "A"]); p2 <- mean(d$cure[d$arm == "B"])
  se <- sqrt(p1 * (1 - p1) / 120 + p2 * (1 - p2) / 110)
  pl <- pnorm((p1 - p2 + 0.1) / se, lower.tail = FALSE)
  pu <- pnorm((p1 - p2 - 0.1) / se)
  expect_equal(r$primary_test$p.value, max(pl, pu))
  expect_equal(c(r$primary_test$conf.low, r$primary_test$conf.high), p1 - p2 + c(-1, 1) * qnorm(0.95) * se)
  expect_match(r$primary_test$method, "proportions")
})

test_that("testequivalence validates the margin", {
  d <- data.frame(y = rnorm(10), g = rep(1:2, 5))
  expect_error(testequivalence(y ~ g, d, lower = 1, upper = -1), "smaller")
  expect_error(testequivalence(y ~ g, d, lower = -Inf, upper = Inf), "infinite")
})
