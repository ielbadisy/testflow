test_that("test_repeated uses repeated-measures ANOVA when normality passes", {
  set.seed(10)
  dat <- tibble::tibble(
    id = seq_len(40),
    t1 = rnorm(40, 10, 2),
    t2 = t1 + rnorm(40, 1, 1),
    t3 = t1 + rnorm(40, 2, 1)
  )
  x <- test_repeated(dat, c(t1, t2, t3), id = id)
  long <- tidyr::pivot_longer(dat, c(t1, t2, t3), names_to = "time", values_to = "value")
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

  long <- tidyr::pivot_longer(cbind(id = factor(seq_len(n)), wide), c(t1, t2, t3), names_to = "time", values_to = "value")
  ref <- summary(stats::aov(value ~ time + Error(id / time), data = long))[["Error: id:time"]][[1]]
  expect_equal(x$primary_test$statistic[1], unname(ref[["F value"]][1]))
  expect_equal(x$primary_test$parameter[1], unname(ref[["Df"]][1]))
  expect_equal(x$primary_test$df.error[1], unname(ref[["Df"]][2]))
  expect_equal(x$primary_test$p.value[1], unname(ref[["Pr(>F)"]][1]))
})
