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
  ref <- summary(stats::aov(value ~ time + Error(id / time), data = long))[["Error: Within"]][[1]][["Pr(>F)"]][1]
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
