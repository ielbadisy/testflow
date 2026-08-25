test_that("test_paired computes differences and returns plot", {
  set.seed(3)
  dat <- tibble::tibble(before = rnorm(50, 10, 2), after = before - rnorm(50, 1, 1))
  x <- test_paired(dat, before, after)
  expect_s3_class(x, "testflow_paired")
  expect_equal(x$data$diff, x$data$after - x$data$before)
  expect_equal(x$recommended$test, "Paired t-test")
  expect_s3_class(plot(x), "ggplot")
  expect_named(x$alternative_tests, c("paired_t_test", "wilcox", "sign"))
  z <- test_paired(after ~ before, data = dat)
  expect_s3_class(z, "testflow_paired")
  expect_equal(z$data$diff, z$data$after - z$data$before)
  expect_true("effect_size" %in% names(z$primary_test))
  expect_false(is.na(z$primary_test$effect_size[1]))
})

test_that("test_paired's Wilcoxon branch reports a Hodges-Lehmann estimate and CI", {
  set.seed(4)
  dat <- tibble::tibble(before = stats::rexp(30, 1), after = before + stats::rexp(30, 1) - 0.5)
  x <- test_paired(dat, before, after)

  wilcox <- x$alternative_tests$wilcox
  expect_false(is.na(wilcox$estimate[1]))
  expect_false(is.na(wilcox$conf.low[1]))
  expect_false(is.na(wilcox$conf.high[1]))
})
