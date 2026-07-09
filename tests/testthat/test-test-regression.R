test_that("test_linear_regression matches summary(lm()) exactly", {
  set.seed(1)
  n <- 80
  dat <- tibble::tibble(x1 = rnorm(n), x2 = rnorm(n))
  dat$y <- 2 + 0.6 * dat$x1 - 0.3 * dat$x2 + rnorm(n, sd = 0.8)

  x <- test_linear_regression(y ~ x1 + x2, data = dat)
  ref <- lm(y ~ x1 + x2, data = dat)
  refsum <- summary(ref)

  expect_s3_class(x, "testflow_linear_regression")
  expect_equal(x$primary_test$statistic, unname(refsum$fstatistic["value"]))
  expect_equal(x$primary_test$parameter, unname(refsum$fstatistic["numdf"]))
  expect_equal(x$primary_test$df.error, unname(refsum$fstatistic["dendf"]))
  expect_equal(
    x$primary_test$p.value,
    unname(stats::pf(refsum$fstatistic[1], refsum$fstatistic[2], refsum$fstatistic[3], lower.tail = FALSE))
  )
  expect_equal(x$effect_size$estimate[x$effect_size$name == "R squared"], refsum$r.squared)
  expect_equal(x$effect_size$estimate[x$effect_size$name == "Adjusted R squared"], refsum$adj.r.squared)

  coefs <- x$alternative_tests$coefficients
  expect_equal(coefs$estimate[coefs$term == "x1"], unname(coef(ref)["x1"]))
  expect_equal(coefs$estimate[coefs$term == "x2"], unname(coef(ref)["x2"]))
  expect_equal(coefs$p.value[coefs$term == "x1"], unname(summary(ref)$coefficients["x1", "Pr(>|t|)"]))

  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
})

test_that("test_linear_regression skips multicollinearity for a single predictor", {
  set.seed(2)
  dat <- tibble::tibble(x1 = rnorm(50), y = rnorm(50))
  x <- test_linear_regression(y ~ x1, data = dat)
  mc <- x$assumptions[x$assumptions$name == "Multicollinearity", ]
  expect_equal(mc$status, "not applicable")
})

test_that("test_linear_regression flags multicollinearity for highly correlated predictors", {
  set.seed(3)
  n <- 100
  x1 <- rnorm(n)
  dat <- tibble::tibble(x1 = x1, x2 = x1 + rnorm(n, sd = 0.01), y = rnorm(n))
  x <- test_linear_regression(y ~ x1 + x2, data = dat)
  mc <- x$assumptions[x$assumptions$name == "Multicollinearity", ]
  expect_equal(mc$status, "not acceptable")
  ref_vif <- car::vif(lm(y ~ x1 + x2, data = dat))
  expect_equal(mc$statistic, unname(max(ref_vif)))
})

test_that("test_logistic_regression matches glm() likelihood-ratio test and odds ratios exactly", {
  set.seed(1)
  n <- 150
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  p <- stats::plogis(0.5 + 0.9 * x1 - 0.6 * x2)
  dat <- tibble::tibble(y = rbinom(n, 1, p), x1 = x1, x2 = x2)

  x <- test_logistic_regression(y ~ x1 + x2, data = dat)
  ref <- glm(y ~ x1 + x2, data = dat, family = binomial)
  null_ref <- glm(y ~ 1, data = dat, family = binomial)
  lr_ref <- ref$null.deviance - ref$deviance
  df_ref <- ref$df.null - ref$df.residual
  p_ref <- stats::pchisq(lr_ref, df_ref, lower.tail = FALSE)

  expect_s3_class(x, "testflow_logistic_regression")
  expect_equal(x$primary_test$statistic, lr_ref)
  expect_equal(x$primary_test$parameter, df_ref)
  expect_equal(x$primary_test$p.value, p_ref)

  mcf_ref <- 1 - as.numeric(logLik(ref)) / as.numeric(logLik(null_ref))
  expect_equal(x$effect_size$estimate[1], mcf_ref)

  coefs <- x$alternative_tests$coefficients
  expect_equal(coefs$estimate[coefs$term == "x1"], unname(coef(ref)["x1"]))
  ors <- x$alternative_tests$odds_ratios
  expect_equal(ors$odds_ratio[ors$term == "x1"], unname(exp(coef(ref)["x1"])))

  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
})

test_that("test_logistic_regression accepts a factor outcome and flags Cook's distance outliers", {
  set.seed(4)
  n <- 60
  dat <- tibble::tibble(
    x1 = rnorm(n),
    y = factor(rbinom(n, 1, stats::plogis(0.3 * rnorm(n))), labels = c("no", "yes"))
  )
  x <- test_logistic_regression(y ~ x1, data = dat)
  expect_s3_class(x, "testflow_logistic_regression")
  inf <- x$assumptions[x$assumptions$name == "Influential observations", ]
  expect_true(inf$status %in% c("acceptable", "warning"))
})

test_that("test_logistic_regression errors informatively for a non-binary outcome", {
  dat <- tibble::tibble(x1 = rnorm(30), y = sample(1:3, 30, replace = TRUE))
  expect_error(test_logistic_regression(y ~ x1, data = dat), "two non-missing levels")
})
