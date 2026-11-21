okcol <- function(v) {
  v <- toupper(unique(stats::na.omit(as.character(v))))
  rgb <- grDevices::col2rgb(v)
  # Greys (equal RGB) or one of the Okabe-Ito accents.
  all(rgb[1, ] == rgb[2, ] & rgb[2, ] == rgb[3, ] | substr(v, 1, 7) %in% c("#0072B2", "#D55E00"))
}
cols <- function(p) unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) c(d$colour, d$fill)))
ys <- function(p) {
  b <- ggplot2::ggplot_build(p)
  unlist(lapply(b$layout$panel_params, function(pp) pp$y$get_labels()))
}

test_that("title, subtitle and caption given as NULL are removed", {
  dat <- makecardiodata(80)
  p <- plot(testtwogroups(sbp_3m ~ sex, data = dat), title = NULL, subtitle = NULL, caption = NULL)
  expect_null(p$labels$title); expect_null(p$labels$subtitle); expect_null(p$labels$caption)
  p2 <- plot(testtwogroups(sbp_3m ~ sex, data = dat))
  expect_false(is.null(p2$labels$subtitle))
})

test_that("plot() errors with a clear message when the test ran with plot = FALSE", {
  dat <- makecardiodata(80)
  expect_error(plot(testtwogroups(sbp_3m ~ sex, data = dat, plot = FALSE)), "plot = FALSE")
})

test_that("plots use greys and the two accents only, with no workflow title", {
  dat <- makecardiodata(120)
  b <- transform(MASS::birthwt, smoke = factor(smoke, 0:1, c("no", "yes")), low = factor(low, 0:1, c("no", "yes")))
  d2 <- data.frame(t1 = rep(c(TRUE, FALSE, TRUE, FALSE), c(30, 10, 5, 55)), t2 = rep(c(TRUE, TRUE, FALSE, FALSE), c(30, 10, 5, 55)),
                   dz = rep(c(TRUE, FALSE), c(45, 55)))
  ps <- list(
    testcategorical(low ~ smoke, data = b),
    testproportion(b, "low", success = "yes", p = 0.3),
    testtrend(low ~ race, b),
    testmultinomial(b, "race", p = c(0.5, 0.25, 0.25)),
    testequivalence(bwt ~ smoke, b, lower = -300, upper = 300),
    testblandaltman(transform(dat, a = sbp_3m, b = sbp_3m + rnorm(120, 1, 5)), "a", "b"),
    testroc(b, "lwt", "low"),
    testroccompare(b, "lwt", "age", "low"),
    testdiagnostic(transform(b, ht = factor(ht, 0:1, c("no", "yes"))), "ht", "low"),
    testpaireddiagnostic(d2, "t1", "t2", "dz"),
    testicc(as.data.frame(matrix(rnorm(90), 30) + rnorm(30)), measures = c("V1", "V2", "V3")),
    testcorrelation(bwt ~ lwt, data = b, method = "pearson", covariates = "age")
  )
  for (x in ps) {
    p <- plot(x)
    expect_s3_class(p, "ggplot"); expect_null(p$labels$title)
    expect_true(okcol(cols(p)), info = class(x)[2])
  }
})

test_that("ratio forests use real term names and per-unit rescaling", {
  b <- transform(MASS::birthwt, smoke = factor(smoke, 0:1, c("No", "Yes")))
  x <- testlogisticregression(low ~ age + smoke, b, units = c(age = 10))
  lab <- ys(plot(x))
  expect_true(all(c("age (per 10)", "smoke: Yes") %in% lab))
  fit <- stats::glm(low ~ age + smoke, binomial, b)
  d <- ggplot2::ggplot_build(plot(x))$data
  xs <- unlist(lapply(d, function(z) z$x))
  expect_true(any(abs(10^xs - exp(10 * stats::coef(fit)[["age"]])) < 1e-8))
  x0 <- testlogisticregression(low ~ age + smoke, b)
  expect_equal(x$primary_test, x0$primary_test)
})

test_that("testcox crude = TRUE stores univariable HRs matching coxph and draws both", {
  skip_if_not_installed("survival")
  lung <- transform(survival::lung, sex = factor(sex, 1:2, c("Male", "Female")))
  x <- testcox(survival::Surv(time, status) ~ age + sex, data = lung, crude = TRUE, units = c(age = 10))
  cr <- x$alternative_tests$crude
  ref <- stats::coef(survival::coxph(survival::Surv(time, status) ~ age, data = lung))
  expect_equal(cr$estimate[grepl("^age", cr$term)], exp(10 * unname(ref)), tolerance = 1e-6)
  expect_true(all(c("Adjusted", "Crude") %in% unlist(lapply(ggplot2::ggplot_build(plot(x))$plot$scales$scales,
                                                                function(s) s$get_breaks()))))
  x0 <- testcox(survival::Surv(time, status) ~ age + sex, data = lung)
  expect_null(x0$alternative_tests$crude)
  expect_equal(x$primary_test, x0$primary_test)
})

test_that("stratified forest has a crude row and flags zero-margin strata", {
  ex <- with(esoph, data.frame(agegp = rep(rep(agegp, 2), c(ncases, ncontrols)),
                               alc = rep(rep(as.integer(alcgp) > 2, 2), c(ncases, ncontrols)),
                               case = rep(c(1, 0), c(sum(ncases), sum(ncontrols)))))
  x <- teststratified(case ~ alc | agegp, ex, event = "1")
  p <- plot(x); lab <- ys(p)
  expect_true(all(c("Crude", "MH pooled") %in% lab))
  d <- ggplot2::ggplot_build(p)$data
  pr <- do.call(rbind, lapply(d, function(z) if ("xmin" %in% names(z)) z[, c("x", "xmin", "xmax")]))
  expect_true(any(abs(10^pr$x - x$alternative_tests$crude_or$odds_ratio) < 1e-8))
  expect_true(any(abs(10^pr$x - x$primary_test$estimate) < 1e-8))
  z <- rbind(ex, data.frame(agegp = "none", alc = c(TRUE, FALSE), case = 0))
  pz <- plot(suppressWarnings(teststratified(case ~ alc | agegp, z, event = "1")))
  tx <- unlist(lapply(ggplot2::ggplot_build(pz)$data, function(z) z$label))
  expect_true(any(grepl("not estimable", tx)))
})

test_that("Wilson and Clopper-Pearson intervals in plots match binom", {
  skip_if_not_installed("binom")
  b <- transform(MASS::birthwt, low = factor(low, 0:1, c("no", "yes")))
  x <- testtrend(low ~ race, b)
  d <- ggplot2::ggplot_build(plot(x))$data
  pr <- do.call(rbind, lapply(d, function(z) if ("ymin" %in% names(z)) z[, c("y", "ymin", "ymax")]))
  ev <- tapply(b$low == "yes", b$race, sum); n <- tapply(b$low, b$race, length)
  ref <- binom::binom.confint(ev, n, methods = "wilson")
  expect_equal(pr$ymin, ref$lower, tolerance = 1e-8); expect_equal(pr$ymax, ref$upper, tolerance = 1e-8)
  # The fitted line is the least-squares trend the Cochran-Armitage test targets.
  fl <- d[[1]]; cf <- stats::coef(stats::lm(I(low == "yes") ~ race, b))
  expect_equal(fl$y, cf[1] + cf[2] * fl$x, tolerance = 1e-8)
})

test_that("generic power curves pass through the target and agree with power.t.test and power.prop.test", {
  for (x in list(samplesizecontinuous(design = "parallel", delta = 3, sd = 8, power = 0.9),
                 samplesizebinary(design = "parallel", p1 = 0.3, p2 = 0.15, power = 0.9),
                 samplesizesurvival(hr = 0.7, power = 0.9))) {
    p <- plot(x, type = "curve"); expect_s3_class(p, "ggplot"); expect_null(p$labels$title)
    cd <- p$data
    expect_equal(stats::approx(cd$n, cd$power, sum(x$plot_data$value[x$plot_data$type == "raw"]))$y, 0.9, tolerance = 0.01)
  }
  x <- samplesizecontinuous(design = "parallel", delta = 3, sd = 8, power = 0.9)
  cd <- plot(x, type = "curve")$data
  ref <- vapply(cd$n / 2, function(m) stats::power.t.test(m, 3, 8)$power, 1)
  expect_equal(cd$power, ref, tolerance = 0.01)
  x <- samplesizebinary(design = "parallel", p1 = 0.3, p2 = 0.15, power = 0.9)
  cd <- plot(x, type = "curve")$data
  ref <- vapply(cd$n / 2, function(m) stats::power.prop.test(m, 0.3, 0.15)$power, 1)
  expect_equal(cd$power, ref, tolerance = 0.02)
})

test_that("partial correlation plot shows the residuals whose correlation is the reported estimate", {
  b <- MASS::birthwt
  x <- testcorrelation(bwt ~ lwt, data = b, method = "pearson", covariates = "age")
  d <- ggplot2::ggplot_build(plot(x))$data[[1]]
  expect_equal(stats::cor(d$x, d$y), x$primary_test$estimate, tolerance = 1e-10)
  ref <- stats::cor(stats::resid(stats::lm(bwt ~ age, b)), stats::resid(stats::lm(lwt ~ age, b)))
  expect_equal(x$primary_test$estimate, ref, tolerance = 1e-10)
})

test_that("testsurvival risk table and testrepeated time labels", {
  skip_if_not_installed("survival")
  lung <- transform(survival::lung, sex = factor(sex, 1:2, c("Male", "Female")))
  x <- testsurvival(survival::Surv(time, status) ~ sex, data = lung, times = c(0, 250, 500, 750))
  tx <- unlist(lapply(ggplot2::ggplot_build(plot(x))$data, function(z) z$label))
  expect_true("At risk" %in% tx)
  ref <- summary(survival::survfit(survival::Surv(time, status) ~ sex, lung), times = c(0, 250, 500, 750))$n.risk
  expect_true(all(as.character(ref) %in% tx))
  x0 <- testsurvival(survival::Surv(time, status) ~ sex, data = lung, risktable = FALSE)
  expect_false("At risk" %in% unlist(lapply(ggplot2::ggplot_build(plot(x0))$data, function(z) z$label)))
  expect_equal(x$primary_test, x0$primary_test)
  w <- data.frame(id = 1:20, m1 = rnorm(20), m2 = rnorm(20), m3 = rnorm(20))
  p <- plot(testrepeated(w, c("m1", "m2", "m3"), id = "id", times = c(0, 6, 12)))
  expect_equal(as.vector(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()), c("0", "6", "12"))
  expect_error(testrepeated(w, c("m1", "m2", "m3"), id = "id", times = 1:2), "times")
})

test_that("paired categorical and repeated categorical plots use real labels", {
  d <- data.frame(b = rep(c("poor", "good", "poor", "good"), c(16, 4, 14, 20)),
                  a = rep(c("poor", "poor", "good", "good"), c(16, 4, 14, 20)))
  p <- plot(testpairedcategorical(d, "b", "a"))
  expect_equal(p$labels$x, "a"); expect_equal(p$labels$y, "b")
  expect_true(all(c("poor", "good") %in% ys(p)))
  r <- data.frame(id = 1:30, t1 = rep(c("no", "yes"), c(20, 10)), t2 = rep(c("no", "yes"), c(12, 18)),
                  t3 = rep(c("no", "yes"), c(8, 22)))
  p <- plot(testrepeatedcategorical(r, c(t1, t2, t3), id = id, event = "yes"))
  expect_equal(p$labels$y, "Percent yes")
  expect_equal(ggplot2::ggplot_build(p)$data[[2]]$y, c(10, 18, 22) / 30)
})
