test_that("test_survival matches survdiff() and a companion coxph() exactly", {
  set.seed(1)
  n <- 100
  dat <- tibble::tibble(
    time = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    group = rep(c("A", "B"), each = n / 2)
  )

  x <- test_survival(survival::Surv(time, status) ~ group, data = dat)
  ref_sd <- survival::survdiff(survival::Surv(time, status) ~ group, data = dat)
  ref_cox <- survival::coxph(survival::Surv(time, status) ~ group, data = dat)
  ref_hr <- unname(exp(coef(ref_cox)))
  ref_ci <- unname(exp(confint(ref_cox)))

  expect_s3_class(x, "testflow_survival")
  expect_equal(x$primary_test$statistic, ref_sd$chisq)
  expect_equal(x$primary_test$parameter, length(ref_sd$n) - 1)
  expect_equal(x$primary_test$p.value, ref_sd$pvalue)
  expect_equal(x$effect_size$estimate[1], ref_hr)
  expect_equal(c(x$primary_test$conf.low, x$primary_test$conf.high), as.numeric(ref_ci))
  expect_equal(x$effect_size$name[1], "Hazard ratio")

  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
  expect_true(grepl("Hazard ratio", report(x), fixed = TRUE))
})

test_that("test_survival requires exactly two groups", {
  set.seed(2)
  n <- 60
  dat <- tibble::tibble(
    time = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    group = sample(c("A", "B", "C"), n, replace = TRUE)
  )
  expect_error(test_survival(survival::Surv(time, status) ~ group, data = dat), "exactly two")
})

test_that("test_cox matches coxph()/cox.zph() exactly for a multi-predictor model", {
  set.seed(1)
  n <- 150
  dat <- tibble::tibble(
    time = rexp(n, 0.08),
    status = rbinom(n, 1, 0.75),
    age = rnorm(n, 60, 10),
    group = rep(c("A", "B"), each = n / 2)
  )

  x <- test_cox(survival::Surv(time, status) ~ age + group, data = dat)
  ref <- survival::coxph(survival::Surv(time, status) ~ age + group, data = dat)
  refsum <- summary(ref)
  ref_zph <- survival::cox.zph(ref)

  expect_s3_class(x, "testflow_cox")
  expect_equal(x$primary_test$statistic, unname(refsum$logtest["test"]))
  expect_equal(x$primary_test$parameter, unname(refsum$logtest["df"]))
  expect_equal(x$primary_test$p.value, unname(refsum$logtest["pvalue"]))
  expect_equal(x$effect_size$estimate[1], unname(refsum$concordance["C"]))
  expect_equal(x$effect_size$name[1], "Concordance index")

  coefs <- x$alternative_tests$hazard_ratios
  expect_equal(coefs$hazard_ratio[coefs$term == "age"], unname(exp(coef(ref)["age"])))

  ph <- x$assumptions[x$assumptions$name == "Proportional hazards", ]
  expect_equal(ph$statistic, unname(ref_zph$table["GLOBAL", "chisq"]))
  expect_equal(ph$p_value, unname(ref_zph$table["GLOBAL", "p"]))

  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
})

test_that("test_cox works with a single predictor", {
  set.seed(3)
  n <- 80
  dat <- tibble::tibble(
    time = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    age = rnorm(n, 60, 10)
  )
  x <- test_cox(survival::Surv(time, status) ~ age, data = dat)
  expect_s3_class(x, "testflow_cox")
  expect_equal(nrow(x$alternative_tests$hazard_ratios), 1)
})

test_that("effect sizes without a magnitude label still print and report the estimate", {
  set.seed(1)
  n <- 100
  dat <- tibble::tibble(
    time = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    group = rep(c("A", "B"), each = n / 2)
  )
  x <- test_survival(survival::Surv(time, status) ~ group, data = dat)
  expect_true(is.na(x$effect_size$magnitude[1]))
  expect_false(is.na(x$effect_size$estimate[1]))
  printed <- capture.output(print(x))
  expect_true(any(grepl("Hazard ratio", printed, fixed = TRUE)))
  expect_false(any(grepl("not reported", printed, fixed = TRUE)))
})
