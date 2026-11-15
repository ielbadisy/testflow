test_that("samplesizeprediction matches pmsampsize for a binary outcome", {
  skip_if_not_installed("pmsampsize")
  ours <- samplesizeprediction("binary", parameters = 24, csrsquared = 0.288, prevalence = 0.174)
  ref <- utils::capture.output(pm <- pmsampsize::pmsampsize(type = "b", csrsquared = 0.288, parameters = 24, prevalence = 0.174))
  expect_equal(ours$n, pm$sample_size)
  expect_equal(ours$criteria$n, unname(pm$results_table[, "Samp_size"]))
  expect_equal(ours$criteria$shrinkage, unname(pm$results_table[, "Shrinkage"]), tolerance = 1e-3)
})

test_that("samplesizeprediction matches pmsampsize with Nagelkerke R-squared", {
  skip_if_not_installed("pmsampsize")
  ours <- samplesizeprediction("binary", parameters = 10, nagrsquared = 0.2, prevalence = 0.3)
  ref <- utils::capture.output(pm <- pmsampsize::pmsampsize(type = "b", nagrsquared = 0.2, parameters = 10, prevalence = 0.3))
  expect_equal(ours$n, pm$sample_size)
})

test_that("samplesizeprediction matches pmsampsize for a survival outcome", {
  skip_if_not_installed("pmsampsize")
  ours <- samplesizeprediction("survival", parameters = 25, csrsquared = 0.051, rate = 0.065, timepoint = 2, meanfup = 2.07)
  ref <- utils::capture.output(pm <- pmsampsize::pmsampsize(type = "s", csrsquared = 0.051, parameters = 25, rate = 0.065, timepoint = 2, meanfup = 2.07))
  expect_equal(ours$n, pm$sample_size)
  expect_equal(ours$criteria$n, unname(pm$results_table[, "Samp_size"]))
})

test_that("samplesizeprediction matches pmsampsize for a continuous outcome", {
  skip_if_not_installed("pmsampsize")
  ours <- samplesizeprediction("continuous", parameters = 25, rsquared = 0.2, intercept = 1.9, sd = 0.6)
  ref <- utils::capture.output(pm <- pmsampsize::pmsampsize(type = "c", rsquared = 0.2, parameters = 25, intercept = 1.9, sd = 0.6))
  expect_equal(ours$n, pm$sample_size)
  expect_equal(ours$criteria$n, unname(pm$results_table[, "Samp_size"]))
})

test_that("samplesizeprediction validates inputs", {
  expect_error(samplesizeprediction("binary", parameters = 10, prevalence = 0.3), "exactly one")
  expect_error(samplesizeprediction("binary", parameters = 10, csrsquared = 0.9, prevalence = 0.1), "maximum")
  expect_s3_class(samplesizeprediction("binary", parameters = 10, csrsquared = 0.1, prevalence = 0.3), "sample_size")
})
