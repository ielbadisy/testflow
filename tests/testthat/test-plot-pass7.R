labs_of <- function(p) ggplot2::ggplot_build(p)$plot$labels
ylabels <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()

test_that("testequivalence takes xlab and units", {
  set.seed(1); d <- data.frame(y = rnorm(40), g = rep(c("F", "M"), 20))
  expect_equal(labs_of(plot(testequivalence(y ~ g, d, lower = -3, upper = 3, units = "kg")))$x,
               "Mean difference (F - M), kg")
  expect_equal(labs_of(plot(testequivalence(y ~ g, d, lower = -3, upper = 3, xlab = "Change")))$x, "Change")
})

test_that("testpaired group facets the plot and leaves the tests unchanged", {
  a <- MASS::anorexia
  x <- testpaired(Postwt ~ Prewt, data = a, group = "Treat"); x0 <- testpaired(Postwt ~ Prewt, data = a)
  expect_s3_class(plot(x)$facet, "FacetWrap"); expect_equal(x$primary_test, x0$primary_test)
  expect_equal(levels(plot(x)$data$Treat), levels(a$Treat))
})

test_that("testcategorical type = 'residuals' plots the adjusted standardised residuals", {
  d <- data.frame(y = rep(c("a", "b", "c"), c(30, 20, 10)), g = rep(c("x", "y"), 30))
  x <- testcategorical(y ~ g, data = d); p <- plot(x, type = "residuals")
  expect_equal(sort(p$data$stdres), sort(as.vector(suppressWarnings(chisq.test(table(d$y, d$g)))$stdres)))
  expect_null(labs_of(plot(x, type = "residuals", caption = NULL))$caption)
  expect_identical(plot(x), x$plot)
})

test_that("testrepeatedcategorical event_label names the plotted level", {
  set.seed(2); w <- data.frame(id = 1:30, t1 = rbinom(30, 1, .3), t2 = rbinom(30, 1, .5), t3 = rbinom(30, 1, .7))
  p <- plot(testrepeatedcategorical(w, t1:t3, id = "id", event = 1, event_label = "good status"))
  expect_equal(labs_of(p)$y, "Percent good status"); expect_match(labs_of(p)$caption, "good status")
})

test_that("testrepeated takes xlab, ylab and legendtitle", {
  set.seed(3); d <- data.frame(id = 1:20, arm = rep(c("A", "B"), 10), m0 = rnorm(20), m1 = rnorm(20), m2 = rnorm(20))
  l <- labs_of(plot(testrepeated(d, m0:m2, id = id, between = arm, xlab = "Month", ylab = "DBP", legendtitle = NULL)))
  expect_equal(c(l$x, l$y), c("Month", "DBP")); expect_null(l$colour)
  expect_equal(labs_of(plot(testrepeated(d, m0:m2, id = id, between = arm)))$colour, "arm")
})

test_that("multinomial forest panels follow the outcome factor order", {
  d <- transform(iris, s = factor(ifelse(Sepal.Length < 5.5, "Small", ifelse(Sepal.Length < 6.5, "Mid", "Big")),
                                  c("Small", "Mid", "Big")))
  p <- plot(testmultinomialregression(s ~ Petal.Width, d))
  expect_equal(levels(p$data$panel), c("Mid", "Big"))
})

test_that("forest labels and unitlabels relabel terms", {
  b <- MASS::birthwt; b$smoke <- factor(b$smoke, 0:1, c("No", "Yes"))
  p <- plot(testlogisticregression(low ~ age + smoke, b, units = c(age = 5), labels = c(smoke = "Smoking"),
                                   unitlabels = c(age = "years")))
  expect_setequal(ylabels(p), c("age (per 5 years)", "Smoking: Yes"))
  p <- plot(testlogisticregression(low ~ age + smoke, b, labels = c(smokeYes = "Smoker", age = "Age")))
  expect_setequal(ylabels(p), c("Age", "Smoker"))
})

test_that("testcox units rescale the HR tables, crude is on the HR scale, test in the caption", {
  skip_if_not_installed("survival")
  lung <- transform(survival::lung, sex = factor(sex, 1:2, c("Male", "Female")))
  f <- survival::Surv(time, status) ~ age + sex
  x <- testcox(f, data = lung, units = c(age = 10), crude = TRUE)
  fit <- survival::coxph(f, data = lung); ci <- confint(fit)
  hr <- x$alternative_tests$hazard_ratios
  expect_equal(hr$hazard_ratio, exp(c(10, 1) * coef(fit)), ignore_attr = TRUE, tolerance = 1e-8)
  expect_equal(hr$conf.low, exp(c(10, 1) * ci[, 1]), ignore_attr = TRUE, tolerance = 1e-8)
  expect_equal(hr$per, c(10, 1))
  cr <- x$alternative_tests$crude
  ref <- survival::coxph(survival::Surv(time, status) ~ sex, data = lung)
  expect_equal(cr$estimate[2], exp(unname(coef(ref))), tolerance = 1e-8)
  expect_equal(cr$conf.high[2], exp(unname(confint(ref)[, 2])), tolerance = 1e-8)
  pd <- plot(x)$data
  expect_equal(pd$est[pd$src == "Adjusted"], hr$hazard_ratio, tolerance = 1e-8)
  l <- labs_of(plot(x, title = NULL)); expect_null(l$subtitle); expect_match(l$caption, "Likelihood-ratio")
  x0 <- testcox(f, data = lung)
  expect_equal(x0$alternative_tests$hazard_ratios$hazard_ratio, exp(coef(fit)), ignore_attr = TRUE)
  expect_false("per" %in% names(x0$alternative_tests$hazard_ratios))
})

test_that("sample size vary draws one curve per value, exact t power", {
  ss <- samplesizecontinuous(design = "parallel", delta = 3, sd = 2.7, power = 0.9, method = "t")
  p <- plot(ss, vary = list(delta = c(2, 3)))
  tg <- p$layers[[3]]$data
  ref <- vapply(c(2, 3), function(d) 2 * ceiling(power.t.test(delta = d, sd = 2.7, power = 0.9)$n), 0)
  expect_equal(tg$n, ref)
  cv <- p$data[p$data$v == levels(p$data$v)[2], ]; i <- which.min(abs(cv$n - 40))
  expect_equal(cv$power[i], power.t.test(n = cv$n[i] / 2, delta = 3, sd = 2.7, strict = TRUE)$power, tolerance = 1e-6)
  pb <- plot(samplesizebinary(p1 = 0.3, p2 = 0.15), vary = list(p2 = c(0.1, 0.15)))
  expect_equal(nlevels(pb$data$v), 2)
  expect_error(plot(ss, vary = c(2, 3)), "named list")
  # Default curve under method = "t" is the exact noncentral t power.
  cd <- plot(ss, type = "curve")$data
  expect_equal(cd$power[5], power.t.test(n = cd$n[5] / 2, delta = 3, sd = 2.7, strict = TRUE)$power, tolerance = 1e-6)
})

test_that("plotclusteradjust matches samplesizeclusteradjust", {
  p <- plotclusteradjust(120, m = c(10, 20), rho = c(0, 0.05))
  expect_equal(p$data$n, c(120, samplesizeclusteradjust(120, 10, 0.05), 120, samplesizeclusteradjust(120, 20, 0.05)))
})
