test_that("test_two_groups returns expected object and Student t-test", {
  set.seed(1)
  dat <- tibble::tibble(y = c(rnorm(40, 0, 1), rnorm(40, 0.5, 1)), g = rep(c("a", "b"), each = 40))
  x <- test_two_groups(dat, y, g)
  expect_s3_class(x, "testflow_two_groups")
  expect_equal(x$recommended$test, "Student independent t-test")
  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
  z <- test_two_groups(y ~ g, data = dat)
  expect_s3_class(z, "testflow_two_groups")
  expect_true("null_hypothesis" %in% names(z$primary_test))
  expect_true("conf.low" %in% names(z$primary_test))
  expect_false(is.na(z$primary_test$conf.low[1]))
  p <- dat |> test_two_groups(y ~ g)
  expect_s3_class(p, "testflow_two_groups")
})

test_that("test_two_groups recommends Welch when variances differ", {
  set.seed(2)
  dat <- tibble::tibble(y = c(rnorm(60, 0, 1), rnorm(60, 0, 4)), g = rep(c("a", "b"), each = 60))
  x <- test_two_groups(dat, y, g)
  expect_equal(x$recommended$test, "Welch t-test")
})

test_that("test_two_groups recommends Wilcoxon when normality fails", {
  dat <- tibble::tibble(y = c(rexp(35), rexp(35, rate = 0.5)), g = rep(c("a", "b"), each = 35))
  x <- test_two_groups(dat, y, g)
  expect_equal(x$recommended$test, "Wilcoxon rank-sum test")
})
