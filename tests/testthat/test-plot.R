test_that("plot.testflow returns stored plot", {
  dat <- make_cardio_data(80)
  x <- test_groups(dat, sbp_3m, treatment)
  expect_s3_class(plot(x), "ggplot")
})

test_that("correlation matrix plot shows cell labels", {
  dat <- make_cardio_data(80)
  x <- suppressWarnings(test_correlation_matrix(dat, c(age, sbp_3m, ldl), plot = TRUE))
  p <- plot(x)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Correlation matrix")
  expect_true(any(vapply(p$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))))
})

test_that("outlier plot shows IQR fences and highlighted points", {
  dat <- make_cardio_data(80)
  x <- suppressWarnings(test_outliers(c(sbp_3m, ldl, crp), data = dat, plot = TRUE))
  p <- plot(x)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Outlier screening")
  expect_true(any(vapply(p$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1))))
  expect_false(is.null(p$facet$params$facets))
})
