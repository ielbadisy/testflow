test_that("one-sample workflow matches base R reference calculations", {
  xdat <- tibble::tibble(y = c(10, 12, 11, 13, 12, 14, 13, 12))
  x <- test_one_sample(xdat, y, mu = 11)
  ref <- stats::t.test(xdat$y, mu = 11)
  expect_equal(unname(x$primary_test$statistic[1]), unname(ref$statistic), tolerance = 1e-12)
  expect_equal(x$primary_test$p.value[1], ref$p.value, tolerance = 1e-12)
  expect_equal(x$primary_test$conf.low[1], ref$conf.int[1], tolerance = 1e-12)
  expect_equal(x$primary_test$conf.high[1], ref$conf.int[2], tolerance = 1e-12)
})

test_that("proportion workflow matches binomial reference", {
  xdat <- tibble::tibble(y = c("yes", "yes", "no", "yes", "no", "yes", "yes", "no"))
  expect_warning(
    x <- test_proportion(xdat, y, success = "yes", p = 0.5),
    "approximation is borderline"
  )
  ref <- stats::binom.test(sum(xdat$y == "yes"), nrow(xdat), p = 0.5)
  expect_equal(x$primary_test$p.value[1], ref$p.value, tolerance = 1e-12)
  expect_equal(x$primary_test$conf.low[1], ref$conf.int[1], tolerance = 1e-12)
  expect_equal(x$primary_test$conf.high[1], ref$conf.int[2], tolerance = 1e-12)
})

test_that("correlation workflows match base cor.test references", {
  dat <- tibble::tibble(x = c(1, 2, 3, 4, 5, 6), y = c(2, 1, 4, 3, 6, 5))
  pear <- test_correlation(y ~ x, data = dat, method = "pearson")
  spear <- test_correlation(y ~ x, data = dat, method = "spearman")
  kend <- test_correlation(y ~ x, data = dat, method = "kendall")
  ref_pear <- stats::cor.test(dat$x, dat$y, method = "pearson")
  ref_spear <- stats::cor.test(dat$x, dat$y, method = "spearman", exact = FALSE)
  ref_kend <- stats::cor.test(dat$x, dat$y, method = "kendall", exact = FALSE)
  expect_equal(pear$primary_test$p.value[1], ref_pear$p.value, tolerance = 1e-12)
  expect_equal(spear$primary_test$p.value[1], ref_spear$p.value, tolerance = 1e-12)
  expect_equal(kend$primary_test$p.value[1], ref_kend$p.value, tolerance = 1e-12)
})

test_that("categorical workflow matches chi-square and Fisher references", {
  dat <- tibble::tibble(x = rep(c("a", "b"), each = 30), y = rep(c("yes", "no"), 30))
  chi <- test_categorical(x ~ y, data = dat)
  ref_chi <- suppressWarnings(stats::chisq.test(table(dat$x, dat$y), correct = FALSE))
  expect_equal(chi$primary_test$p.value[1], ref_chi$p.value, tolerance = 1e-12)
  small <- tibble::tibble(x = c("a", "a", "a", "b"), y = c("yes", "yes", "no", "no"))
  fish <- test_categorical(x ~ y, data = small)
  ref_fisher <- stats::fisher.test(table(small$x, small$y))
  expect_equal(fish$primary_test$p.value[1], ref_fisher$p.value, tolerance = 1e-12)
})

test_that("factorial workflow matches base aov output for a main effect", {
  dat <- tibble::tibble(
    y = c(10, 11, 12, 13, 14, 15, 16, 17),
    a = factor(rep(c("low", "high"), each = 4)),
    b = factor(rep(c("x", "y"), times = 4))
  )
  x <- test_factorial(y ~ a * b, data = dat)
  ref <- broom::tidy(stats::aov(y ~ a * b, data = dat))
  expect_equal(x$primary_test$p.value[1], ref$p.value[1], tolerance = 1e-12)
  expect_equal(x$primary_test$statistic[1], ref$statistic[1], tolerance = 1e-12)
})
