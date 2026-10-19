test_that("test_repeated uses repeated-measures ANOVA when normality passes", {
  set.seed(10)
  dat <- tibble::tibble(
    id = seq_len(40),
    t1 = rnorm(40, 10, 2),
    t2 = t1 + rnorm(40, 1, 1),
    t3 = t1 + rnorm(40, 2, 1)
  )
  x <- test_repeated(dat, c(t1, t2, t3), id = id)
  long <- data.frame(id = rep(dat$id, 3), time = rep(c("t1", "t2", "t3"), each = nrow(dat)), value = c(dat$t1, dat$t2, dat$t3))
  long$id <- factor(long$id)
  ref <- summary(stats::aov(value ~ time + Error(id / time), data = long))[["Error: id:time"]][[1]][["Pr(>F)"]][1]
  expect_s3_class(x, "testflow_repeated")
  expect_equal(x$recommended$test, "Repeated-measures ANOVA")
  expect_equal(x$primary_test$p.value[1], ref, tolerance = 1e-10)
  expect_equal(x$effect_size$name[1], "eta squared")
  expect_true(inherits(x$posthoc, "data.frame"))
})

test_that("test_repeated uses Friedman and paired Wilcoxon when normality fails", {
  dat <- tibble::tibble(
    id = seq_len(30),
    t1 = rexp(30, 1),
    t2 = rexp(30, 0.8),
    t3 = rexp(30, 0.6)
  )
  x <- test_repeated(dat, c(t1, t2, t3), id = id)
  ref <- stats::friedman.test(as.matrix(dat[, c("t1", "t2", "t3")]))
  expect_equal(x$recommended$test, "Friedman test")
  expect_equal(x$primary_test$p.value[1], ref$p.value)
  expect_equal(unique(x$posthoc$method), "paired Wilcoxon signed-rank test")
})

test_that("test_repeated coerces id to a factor so Error() builds the correct stratum", {
  set.seed(42)
  n <- 15
  wide <- data.frame(t1 = rnorm(n), t2 = rnorm(n), t3 = rnorm(n))
  x <- test_repeated(wide, measures = c(t1, t2, t3))
  expect_s3_class(x$data[[x$id]], "factor")

  long <- data.frame(id = factor(rep(seq_len(n), 3)), time = rep(c("t1", "t2", "t3"), each = n), value = c(wide$t1, wide$t2, wide$t3))
  ref <- summary(stats::aov(value ~ time + Error(id / time), data = long))[["Error: id:time"]][[1]]
  expect_equal(x$primary_test$statistic[1], unname(ref[["F value"]][1]))
  expect_equal(x$primary_test$parameter[1], unname(ref[["Df"]][1]))
  expect_equal(x$primary_test$df.error[1], unname(ref[["Df"]][2]))
  expect_equal(x$primary_test$p.value[1], unname(ref[["Pr(>F)"]][1]))
})
