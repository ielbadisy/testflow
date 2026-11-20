set.seed(3)
bdat <- data.frame(y = c(rexp(25), rexp(30, .5)), g = rep(c("a", "b"), c(25, 30)))

test_that("testbootstrap intervals match boot.ci given the same replicates", {
  skip_if_not_installed("boot")
  r <- testbootstrap(y ~ g, bdat, statistic = "median", B = 999, seed = 11, plot = FALSE)
  f <- function(d, i) { x <- d[i, ]; median(x$y[x$g == "b"]) - median(x$y[x$g == "a"]) }
  set.seed(11); b <- boot::boot(bdat, f, R = 999, strata = factor(bdat$g)); b$t[, 1] <- r$replicates
  ci <- boot::boot.ci(b, type = c("norm", "perc", "basic", "bca"), L = boot::empinf(b, type = "jack"))
  iv <- r$alternative_tests$intervals
  ref <- rbind(ci$normal[2:3], ci$percent[4:5], ci$basic[4:5], ci$bca[4:5])
  expect_equal(unname(as.matrix(iv[match(c("normal", "percentile", "basic", "bca"), iv$type), c("conf.low", "conf.high")])),
               ref, tolerance = 1e-8)
  expect_equal(r$primary_test$estimate, f(bdat, seq_len(nrow(bdat))))
  expect_true(all(c("estimate", "std.error", "bias", "conf.low", "conf.high", "B") %in% names(r$primary_test)))
  r1 <- testbootstrap(y ~ 1, bdat, B = 500, seed = 2, plot = FALSE)
  b1 <- boot::boot(bdat$y, function(x, i) mean(x[i]), R = 500); b1$t[, 1] <- r1$replicates
  c1 <- boot::boot.ci(b1, type = "bca", L = boot::empinf(b1, type = "jack"))
  expect_equal(r1$primary_test$conf.low, c1$bca[4], tolerance = 1e-8)
  expect_output(print(r1))
})

test_that("testpermutation exact p matches coin", {
  skip_if_not_installed("coin")
  set.seed(1); d <- data.frame(y = c(rnorm(7), rnorm(6, 1)), g = factor(rep(c("a", "b"), c(7, 6))))
  p <- testpermutation(y ~ g, d, plot = FALSE)
  expect_true(p$alternative_tests$permutation$exact)
  expect_equal(p$primary_test$p.value, coin::pvalue(coin::oneway_test(y ~ g, d, distribution = "exact")), tolerance = 1e-8)
  mc <- testpermutation(y ~ g, bdat, B = 999, seed = 1, plot = FALSE)
  expect_false(mc$alternative_tests$permutation$exact)
  expect_equal(mc$primary_test$p.value, (1 + mc$alternative_tests$permutation$n_extreme) / 1000)
})

test_that("testonesample median CI matches DescTools exact", {
  skip_if_not_installed("DescTools")
  set.seed(4); x <- rnorm(30)
  m <- testonesample(data.frame(x = x), x, plot = FALSE)$alternative_tests$median_ci
  ref <- DescTools::MedianCI(x, method = "exact")
  expect_equal(c(m$estimate, m$conf.low, m$conf.high), unname(ref[1:3]), tolerance = 1e-10)
})

test_that("testproportion exact = 'central' doubles the smaller tail (4/55 vs 0.27)", {
  d <- data.frame(y = c(rep("y", 4), rep("n", 51)))
  a <- testproportion(d, y, "y", p = .27, plot = FALSE)
  b <- testproportion(d, y, "y", p = .27, plot = FALSE, exact = "central")
  expect_equal(a$alternative_tests$exact_binomial$p.value, binom.test(4, 55, .27)$p.value)
  expect_equal(b$alternative_tests$exact_binomial$p.value, 2 * pbinom(4, 55, .27))
})

test_that("testcategorical 2x2 effect measures and residuals (indomethacin repro)", {
  skip_if_not_installed("DescTools")
  tab <- matrix(c(255, 52, 268, 27), 2, dimnames = list(outcome = c("0_no", "1_yes"), rx = c("0_placebo", "1_indomethacin")))
  d <- as.data.frame(as.table(tab)); d <- d[rep(1:4, d$Freq), 1:2]
  r <- testcategorical(outcome ~ rx, d, plot = FALSE)
  em <- r$alternative_tests$effect_measures
  expect_equal(em$measure, c("RD", "RR", "OR", "NNT"))
  rd <- DescTools::BinomDiffCI(27, 295, 52, 307, method = "score")
  expect_equal(unlist(em[1, c("estimate", "conf.low", "conf.high")]), unname(rd[1, ]), tolerance = 1e-6, ignore_attr = TRUE)
  or <- DescTools::OddsRatio(matrix(c(27, 268, 52, 255), 2, byrow = TRUE), conf.level = .95)
  expect_equal(unlist(em[3, c("estimate", "conf.low", "conf.high")]), unname(or), tolerance = 1e-6, ignore_attr = TRUE)
  rr <- 27 / 295 / (52 / 307); se <- sqrt(1 / 27 - 1 / 295 + 1 / 52 - 1 / 307)
  expect_equal(unlist(em[2, c("estimate", "conf.low", "conf.high")]), rr * exp(c(0, -1, 1) * qnorm(.975) * se),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(em$estimate[4], 1 / abs(em$estimate[1]))
  chi <- chisq.test(table(d$outcome, d$rx), correct = FALSE)
  expect_equal(sort(abs(r$posthoc$stdres)), sort(abs(as.vector(chi$stdres))), tolerance = 1e-8)
})

test_that("testtrend matches prop.trend.test and coin::lbl_test", {
  skip_if_not_installed("coin")
  x <- c(10, 15, 22, 30); n <- rep(50, 4)
  d <- data.frame(y = unlist(lapply(1:4, function(i) rep(1:0, c(x[i], n[i] - x[i])))), dose = rep(1:4, n))
  r <- testtrend(y ~ dose, d, plot = FALSE)
  expect_equal(r$primary_test$statistic, unname(prop.trend.test(x, n)$statistic), tolerance = 1e-8)
  tt <- r$alternative_tests$trend_table
  expect_equal(tt$statistic[3], unname(prop.test(x, n)$statistic), tolerance = 1e-8)
  set.seed(1)
  d2 <- data.frame(y = factor(sample(c("lo", "mid", "hi"), 120, TRUE), c("lo", "mid", "hi")), g = factor(sample(letters[1:4], 120, TRUE)))
  r2 <- testtrend(y ~ g, d2, plot = FALSE)
  expect_equal(r2$primary_test$statistic, unname(coin::statistic(coin::lbl_test(table(d2$g, d2$y)))^2), tolerance = 1e-8)
})

test_that("teststratified matches mantelhaen.test and Breslow-Day", {
  skip_if_not_installed("DescTools")
  u <- as.data.frame(UCBAdmissions); u <- u[rep(seq_len(nrow(u)), u$Freq), 1:3]
  s <- teststratified(Admit ~ Gender | Dept, u, plot = FALSE)
  arr <- table(factor(u$Gender, c("Female", "Male")), factor(u$Admit, c("Rejected", "Admitted")), u$Dept)
  mh <- mantelhaen.test(arr)
  expect_equal(s$primary_test$statistic, unname(mh$statistic), tolerance = 1e-8)
  expect_equal(c(s$primary_test$conf.low, s$primary_test$conf.high), as.vector(mh$conf.int), tolerance = 1e-6)
  expect_equal(s$alternative_tests$breslow_day$statistic, unname(DescTools::BreslowDayTest(arr)$statistic), tolerance = 1e-6)
  expect_equal(nrow(s$alternative_tests$stratum_or), 6)
})

test_that("partial correlation matches lm", {
  r <- testcorrelation(mpg ~ hp, mtcars, method = "pearson", covariates = c("wt", "cyl"), plot = FALSE)
  fit <- summary(lm(mpg ~ hp + wt + cyl, mtcars))$coef["hp", ]
  expect_equal(r$primary_test$estimate, cor(resid(lm(mpg ~ wt + cyl, mtcars)), resid(lm(hp ~ wt + cyl, mtcars))))
  expect_equal(c(r$primary_test$statistic, r$primary_test$p.value), unname(fit[3:4]), tolerance = 1e-8)
  expect_false(is.null(r$alternative_tests$marginal))
  m <- suppressWarnings(testcorrelationmatrix(mtcars, c(mpg, hp, disp), method = "pearson", covariates = "wt", plot = FALSE))
  expect_s3_class(m, "testflow")
})

test_that("testroc bootstrap CI for the Youden threshold", {
  set.seed(2); d <- data.frame(m = c(rnorm(60, 1), rnorm(80)), dis = rep(1:0, c(60, 80)))
  a <- testroc(d, m, dis, plot = FALSE)
  b <- testroc(d, m, dis, plot = FALSE, B = 200, seed = 1)
  ot <- b$alternative_tests$optimal_threshold
  expect_null(a$alternative_tests$optimal_threshold$threshold.conf.low)
  expect_true(ot$threshold.conf.low <= ot$threshold & ot$threshold <= ot$threshold.conf.high)
  expect_equal(a$primary_test, b$primary_test)
})

test_that("testpaireddiagnostic matches mcnemar.test", {
  set.seed(1); d <- data.frame(dis = rbinom(200, 1, 0.4))
  d$t1 <- ifelse(d$dis == 1, rbinom(200, 1, 0.8), rbinom(200, 1, 0.1))
  d$t2 <- ifelse(d$dis == 1, rbinom(200, 1, 0.9), rbinom(200, 1, 0.2))
  r <- testpaireddiagnostic(d, t1, t2, dis, plot = FALSE)
  D <- d[d$dis == 1, ]; N <- d[d$dis == 0, ]
  expect_equal(r$primary_test$p.value, c(mcnemar.test(table(D$t1, D$t2))$p.value, mcnemar.test(table(N$t1, N$t2))$p.value))
  expect_equal(r$primary_test$estimate[1], mean(D$t2) - mean(D$t1))
})

test_that("samplesizecorrelation and ANCOVA option", {
  skip_if_not_installed("pwr")
  s <- samplesizecorrelation(.3)
  expect_equal(s$n, ((qnorm(.975) + qnorm(.9)) / atanh(.3))^2 + 3)
  expect_lt(abs(s$n - pwr::pwr.r.test(r = .3, power = .9)$n), 1.5)
  a <- samplesizecontinuous(delta = 5, sd = 10, baseline_cor = .5)
  expect_equal(a$n, samplesizecontinuous(delta = 5, sd = 10 * sqrt(.75))$n)
  expect_equal(samplesizecontinuous(delta = 5, sd = 10)$n, samplesizecontinuous(delta = 5, sd = 10, baseline_cor = 0)$n)
})

test_that("sumtab SMD matches tableone", {
  skip_if_not_installed("tableone")
  d <- transform(mtcars, cyl = factor(cyl), am = factor(am), gear = as.character(gear))
  s <- sumtab(~ mpg + cyl + gear + vs | am, d, smd = TRUE)
  ref <- tableone::ExtractSmd(tableone::CreateTableOne(c("mpg", "cyl", "gear", "vs"), "am", d))
  expect_equal(as.numeric(s$SMD[!is.na(s$SMD)]), round(unname(ref[, 1]), 3))
  expect_error(sumtab(~ mpg | cyl, d, smd = TRUE))
})

test_that("testcox gives per-term PH rows matching cox.zph", {
  lung <- na.omit(survival::lung[, c("time", "status", "age", "sex")])
  x <- testcox(survival::Surv(time, status) ~ age + sex, data = lung, plot = FALSE)
  z <- survival::cox.zph(survival::coxph(survival::Surv(time, status) ~ age + sex, lung))
  ph <- x$assumptions[grepl("^Proportional hazards \\(", x$assumptions$name), ]
  expect_equal(ph$p_value, unname(z$table[c("age", "sex"), "p"]), tolerance = 1e-8)
})

test_that("testlinearregression pretests = FALSE skips residual tests", {
  x <- testlinearregression(mpg ~ wt + hp, mtcars, plot = FALSE, pretests = FALSE)
  a <- x$assumptions[x$assumptions$name %in% c("Normality of residuals", "Homoscedasticity"), ]
  expect_equal(a$status, rep("not checked", 2))
  expect_equal(x$primary_test, testlinearregression(mpg ~ wt + hp, mtcars, plot = FALSE)$primary_test)
})

test_that("testordinalregression matches polr and brant per-term tests", {
  skip_if_not_installed("MASS")
  h <- MASS::housing[rep(1:72, MASS::housing$Freq), ]
  r <- testordinalregression(Sat ~ Infl + Type + Cont, h, plot = FALSE)
  fit <- MASS::polr(Sat ~ Infl + Type + Cont, h, Hess = TRUE)
  or <- r$alternative_tests$odds_ratios
  expect_equal(or$odds_ratio, unname(exp(coef(fit))), tolerance = 1e-6)
  expect_equal(or$conf.low, unname(exp(confint.default(fit)[, 1])), tolerance = 1e-6)
  # brant::brant(fit) per-term X2 (brant 0.3-0)
  expect_equal(r$alternative_tests$proportional_odds$statistic[-1],
               c(0.31733851352, 0.19650458818, 0.46815212458, 2.80475388419, 0.03122600473, 1.28594116758), tolerance = 1e-6)
  expect_output(print(r))
})

test_that("testmultinomialregression matches multinom", {
  skip_if_not_installed("nnet")
  set.seed(5); d <- data.frame(x = rnorm(300), z = rnorm(300))
  d$y <- factor(apply(cbind(0, .5 * d$x, -.5 * d$z + .3) + matrix(-log(-log(runif(900))), 300), 1, which.max))
  r <- testmultinomialregression(y ~ x + z, d, plot = FALSE)
  fit <- nnet::multinom(y ~ x + z, d, trace = FALSE, maxit = 1000)
  co <- r$alternative_tests$coefficients
  expect_equal(co$estimate, as.vector(t(coef(fit))), tolerance = 1e-4)
  expect_equal(co$std.error, as.vector(t(summary(fit)$standard.errors)), tolerance = 1e-3)
  expect_equal(r$alternative_tests$relative_risk_ratios$rrr, exp(co$estimate))
})
