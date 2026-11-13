test_that("plot.testflow returns stored plot", {
  dat <- makecardiodata(80)
  x <- testgroups(dat, sbp_3m, treatment)
  expect_s3_class(plot(x), "ggplot")
  expect_equal(plot(x)$labels$title, x$plot$labels$title)
})

test_that("plot.testflow allows label overrides", {
  dat <- makecardiodata(80)
  x <- testtwogroups(sbp_3m ~ sex, data = dat)
  p <- plot(x, title = "Custom title", subtitle = "Custom subtitle", caption = "Custom caption")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Custom title")
  expect_equal(p$labels$subtitle, "Custom subtitle")
  expect_equal(p$labels$caption, "Custom caption")
  expect_equal(plot(x)$labels$title, x$plot$labels$title)
})

test_that("correlation matrix plot shows cell labels", {
  dat <- makecardiodata(80)
  x <- suppressWarnings(testcorrelationmatrix(dat, c(age, sbp_3m, ldl), plot = TRUE))
  p <- plot(x)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Correlation matrix")
  expect_true(any(vapply(p$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))))
})

test_that("outlier plot shows IQR fences and highlighted points", {
  dat <- makecardiodata(80)
  x <- suppressWarnings(testoutliers(c(sbp_3m, ldl, crp), data = dat, plot = TRUE))
  p <- plot(x)

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Outlier screening")
  expect_true(any(vapply(p$layers, function(layer) inherits(layer$geom, "GeomHline"), logical(1))))
  expect_false(is.null(p$facet$params$facets))
})
