test_that("vary hr recomputes survival_b under proportional hazards, matches separate calls", {
  S <- 0.7; f <- function(h, ...) samplesizesurvival(hr = h, survival_a = S, power = 0.9, accrual_duration = 2,
                                                     follow_up = 3, ...)
  x <- f(0.75, survival_b = S^0.75); hs <- c(0.7, 0.8)
  tg <- plot(x, vary = list(hr = hs))$layers[[3]]$data
  expect_equal(tg$n, vapply(hs, function(h) f(h, survival_b = S^h)$n, 0))
  expect_equal(f(0.75)$n, x$n)
  te <- plot(x, vary = list(hr = hs), what = "events")
  expect_equal(te$layers[[3]]$data$n, vapply(hs, function(h) f(h)$required_events, 0))
  expect_equal(te$labels$x, "Total events")
  expect_error(plot(samplesizebinary(p1 = 0.3, p2 = 0.15), vary = list(p2 = 0.1), what = "events"), "survival")
  # Without survival_a the size is the number of events, as before.
  e <- samplesizesurvival(hr = 0.7, power = 0.9); expect_equal(e$n, e$required_events, tolerance = 1)
})

test_that("vary legend is rounded and takes varylabels", {
  x <- samplesizebinary(p1 = 0.169, p2 = 0.09)
  expect_equal(levels(plot(x, vary = list(p2 = c(0.09152542, 0.12)))$data$v), c("0.0915", "0.12"))
  expect_equal(levels(plot(x, vary = list(p2 = c(0.09, 0.12)), varylabels = c("A", "B"))$data$v), c("A", "B"))
  expect_equal(levels(plot(x, vary = list(p2 = c(0.09, 0.12)), varylabels = c("0.12" = "B"))$data$v), c("0.09", "B"))
})

test_that("custom label on a numeric term keeps the per-unit suffix", {
  b <- MASS::birthwt; yl <- function(...) ggplot2::ggplot_build(plot(testlogisticregression(low ~ age + lwt, b,
    units = c(age = 5), unitlabels = c(age = "years"), ...)))$layout$panel_params[[1]]$y$get_labels()
  expect_setequal(yl(labels = c(age = "Age")), c("Age (per 5 years)", "lwt"))
  expect_setequal(yl(labels = c(age = "Age per 5 years")), c("Age per 5 years", "lwt"))
  expect_setequal(yl(), c("age (per 5 years)", "lwt"))
})
