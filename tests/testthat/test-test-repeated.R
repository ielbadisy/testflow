test_that("testrepeated uses repeated-measures ANOVA when normality passes", {
  set.seed(10)
  dat <- .bt(
    id = seq_len(40),
    t1 = rnorm(40, 10, 2),
    t2 = t1 + rnorm(40, 1, 1),
    t3 = t1 + rnorm(40, 2, 1)
  )
  x <- testrepeated(dat, c(t1, t2, t3), id = id)
  long <- data.frame(id = rep(dat$id, 3), time = rep(c("t1", "t2", "t3"), each = nrow(dat)), value = c(dat$t1, dat$t2, dat$t3))
  long$id <- factor(long$id)
  ref <- summary(stats::aov(value ~ time + Error(id / time), data = long))[["Error: id:time"]][[1]][["Pr(>F)"]][1]
  expect_s3_class(x, "testflow_repeated")
  expect_equal(x$recommended$test, "Repeated-measures ANOVA")
  expect_equal(x$primary_test$p.value[1], ref, tolerance = 1e-10)
  expect_equal(x$effect_size$name[1], "eta squared")
  expect_true(inherits(x$posthoc, "data.frame"))
})

test_that("testrepeated uses Friedman and paired Wilcoxon when normality fails", {
  dat <- .bt(
    id = seq_len(30),
    t1 = rexp(30, 1),
    t2 = rexp(30, 0.8),
    t3 = rexp(30, 0.6)
  )
  x <- testrepeated(dat, c(t1, t2, t3), id = id)
  ref <- stats::friedman.test(as.matrix(dat[, c("t1", "t2", "t3")]))
  expect_equal(x$recommended$test, "Friedman test")
  expect_equal(x$primary_test$p.value[1], ref$p.value)
  expect_equal(unique(x$posthoc$method), "paired Wilcoxon signed-rank test")
})

test_that("testrepeated coerces id to a factor so Error() builds the correct stratum", {
  set.seed(42)
  n <- 15
  wide <- data.frame(t1 = rnorm(n), t2 = rnorm(n), t3 = rnorm(n))
  x <- testrepeated(wide, measures = c(t1, t2, t3))
  expect_s3_class(x$data[[x$id]], "factor")

  long <- data.frame(id = factor(rep(seq_len(n), 3)), time = rep(c("t1", "t2", "t3"), each = n), value = c(wide$t1, wide$t2, wide$t3))
  ref <- summary(stats::aov(value ~ time + Error(id / time), data = long))[["Error: id:time"]][[1]]
  expect_equal(x$primary_test$statistic[1], unname(ref[["F value"]][1]))
  expect_equal(x$primary_test$parameter[1], unname(ref[["Df"]][1]))
  expect_equal(x$primary_test$df.error[1], unname(ref[["Df"]][2]))
  expect_equal(x$primary_test$p.value[1], unname(ref[["Pr(>F)"]][1]))
})

test_that("between factor fits a mixed (split-plot) ANOVA in both repeated workflows", {
  set.seed(11)
  n <- 20
  wide <- data.frame(
    id = seq_len(n),
    grp = factor(rep(c("ctl", "trt"), each = n / 2)),
    t1 = rnorm(n, 10), t2 = rnorm(n, 11), t3 = rnorm(n, 12)
  )
  wide$t3 <- wide$t3 + ifelse(wide$grp == "trt", 2, 0)
  long <- data.frame(
    id = factor(rep(wide$id, 3)), grp = rep(wide$grp, 3),
    time = factor(rep(c("t1", "t2", "t3"), each = n)),
    y = c(wide$t1, wide$t2, wide$t3)
  )
  ref <- summary(aov(y ~ time * grp + Error(id / time), data = long))
  between_tab <- ref[["Error: id"]][[1]]
  within_tab <- ref[["Error: id:time"]][[1]]

  x <- testrepeatedlong(long, y, time, id, between = grp, plot = FALSE)
  m <- x$alternative_tests$mixed_anova
  expect_false(is.null(m))
  expect_setequal(m$term, c("grp", "time", "time:grp"))
  expect_equal(m$statistic[m$term == "grp"], between_tab[["F value"]][1])
  expect_equal(m$p.value[m$term == "time:grp"], within_tab[["Pr(>F)"]][2])
  expect_equal(m$den.df[m$term == "time:grp"], within_tab[["Df"]][3])
  expect_match(report(x), "time:grp", fixed = TRUE)

  w <- testrepeated(wide, c(t1, t2, t3), id = id, between = grp, plot = FALSE)
  mw <- w$alternative_tests$mixed_anova
  expect_equal(mw$statistic[mw$term == "time:grp"], within_tab[["F value"]][2])

  z <- testrepeatedlong(long, y, time, id, plot = FALSE)
  expect_null(z$alternative_tests$mixed_anova)
})
