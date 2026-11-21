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

# Pass 4 plot fixes: sober palette, no workflow titles, labelled axes, means with CI.
geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], "")
builds <- function(p) { b <- ggplot2::ggplot_build(p); expect_true(length(b$data) > 0); b }
okcols <- c("grey15", "grey45", "grey75", "#0072B2", "#D55E00")

test_that("group plots show points, box, mean and CI, no workflow title, effect CI", {
  skip_if_not_installed("MASS")
  b <- transform(MASS::birthwt, smoke = factor(smoke, 0:1, c("no", "yes")))
  p <- plot(testtwogroups(bwt ~ smoke, data = b))
  builds(p)
  expect_null(p$labels$title)
  expect_true(all(c("GeomBoxplot", "GeomPoint", "GeomPointrange") %in% geoms(p)))
  expect_equal(c(p$labels$x, p$labels$y), c("smoke", "bwt"))
  expect_match(p$labels$caption, "Cohen's d = 0.40 \\(95% CI 0.10 to 0.69\\)")
  expect_false(grepl("small|moderate|large", p$labels$caption))
  g <- plot(testgroups(Postwt ~ Treat, data = MASS::anorexia))
  bd <- builds(g)
  expect_true(all(unlist(lapply(bd$data, `[[`, "colour")) %in% okcols))
  expect_false(grepl("workflow", paste(g$labels$title, g$labels$subtitle)))
})

test_that("paired plot puts the first (before) measure first with subject lines", {
  skip_if_not_installed("MASS")
  p <- plot(testpaired(Postwt ~ Prewt, data = subset(MASS::anorexia, Treat == "CBT")))
  b <- builds(p)
  expect_equal(levels(p$data$time), c("Prewt", "Postwt"))
  expect_equal(b$layout$panel_params[[1]]$x$get_labels(), c("Prewt", "Postwt"))
  expect_true(all(c("GeomLine", "GeomPointrange") %in% geoms(p)))
  expect_equal(c(p$labels$x, p$labels$y), c("Measurement", "Value"))
})

test_that("factorial plot shows dodged means with CI, at most three colours for three levels", {
  d <- makecardiodata(80)
  p <- plot(testfactorial(sbp_3m ~ sex * treatment, data = d))
  b <- builds(p)
  expect_null(p$labels$title)
  expect_true("GeomPointrange" %in% geoms(p))
  expect_match(p$labels$subtitle, "sex:treatment p")
  expect_true(all(b$data[[3]]$colour %in% okcols))
})

test_that("repeated plot orders measures and shows between arms", {
  d <- makecardiodata(80)
  x <- testrepeated(d, measures = c(sbp_baseline, sbp_3m, sbp_6m), id = id, between = treatment)
  p <- plot(x)
  b <- builds(p)
  expect_equal(b$layout$panel_params[[1]]$x$get_labels(), c("sbp_baseline", "sbp_3m", "sbp_6m"))
  pr <- b$data[[which(geoms(p) == "GeomPointrange")]]
  expect_equal(nrow(pr), 9)
  expect_equal(length(unique(pr$colour)), 3)
  expect_equal(p$labels$colour, "treatment")
  expect_match(p$labels$subtitle, "Mixed ANOVA")
  q <- plot(testrepeated(d, measures = c(sbp_baseline, sbp_3m, sbp_6m), id = id))
  expect_equal(nrow(builds(q)$data[[which(geoms(q) == "GeomPointrange")]]), 3)
})

test_that("correlation plot gives r with CI in one subtitle and no caption overlap", {
  d <- makecardiodata(80)
  x <- testcorrelation(sbp_3m ~ age, data = d, method = "pearson")
  p <- plot(x)
  builds(p)
  ci <- stats::cor.test(d$age, d$sbp_3m)$conf.int
  expect_match(p$labels$subtitle, paste0("Pearson r = .*95% CI ", format(round(ci[1], 2), nsmall = 2)))
  expect_null(p$labels$title)
  expect_false(grepl("negligible|small", p$labels$subtitle))
})

test_that("plot changes leave the numbers unchanged", {
  skip_if_not_installed("MASS")
  a <- subset(MASS::anorexia, Treat == "CBT")
  expect_equal(testpaired(Postwt ~ Prewt, data = a)$primary_test,
               testpaired(Postwt ~ Prewt, data = a, plot = FALSE)$primary_test)
})
