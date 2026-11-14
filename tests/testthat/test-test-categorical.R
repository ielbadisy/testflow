test_that("testcategorical builds table and computes Cramer's V", {
  dat <- .bt(x = rep(c("a", "b"), each = 30), y = rep(c("yes", "no"), 30))
  z <- testcategorical(dat, x, y)
  expect_s3_class(z, "testflow_categorical")
  expect_equal(z$recommended$test, "Chi-square test of independence")
  expect_equal(z$effect_size$name, "Cramer's V")
  f <- testcategorical(x ~ y, data = dat)
  expect_s3_class(f, "testflow_categorical")
})

test_that("testcategorical uses Fisher with small expected counts", {
  dat <- .bt(x = c("a", "a", "a", "b"), y = c("yes", "yes", "no", "no"))
  z <- testcategorical(dat, x, y)
  expect_equal(z$recommended$test, "Fisher exact test")
})

test_that("testmultinomial converts table counts to numeric internally", {
  dat <- .bt(x = c("a", "a", "b", "c", "c", "c"))
  z <- testmultinomial(dat, x, p = c(1 / 3, 1 / 3, 1 / 3))
  expect_s3_class(z, "testflow_multinomial")
  expect_equal(z$posthoc$expected, rep(2, 3))
})

test_that("testmultinomial adjusts pairwise post-hoc p-values, matching the rest of the package's BH convention", {
  dat <- .bt(x = sample(c("a", "b", "c", "d"), 200, replace = TRUE, prob = c(0.55, 0.2, 0.15, 0.1)))
  z <- suppressWarnings(testmultinomial(dat, x))
  expect_true(all(c("p.adj", "p.adjust.method") %in% names(z$posthoc)))
  expect_equal(z$posthoc$p.adjust.method, rep("BH", nrow(z$posthoc)))
  expect_equal(z$posthoc$p.adj, stats::p.adjust(z$posthoc$p, method = "BH"))
})

test_that("testcategorical offers the Bergsma bias-corrected Cramer's V", {
  dat <- data.frame(cyl = factor(mtcars$cyl), gear = factor(mtcars$gear))
  raw <- testcategorical(cyl ~ gear, data = dat, plot = FALSE)
  adj <- testcategorical(cyl ~ gear, data = dat, plot = FALSE, biascorrect = TRUE)
  tab <- table(dat$cyl, dat$gear)
  n <- sum(tab); r <- nrow(tab); k <- ncol(tab)
  phi2 <- suppressWarnings(chisq.test(tab, correct = FALSE))$statistic / n
  phi2c <- max(0, phi2 - (k - 1) * (r - 1) / (n - 1))
  rc <- r - (r - 1)^2 / (n - 1); kc <- k - (k - 1)^2 / (n - 1)
  expect_equal(adj$effect_size$estimate, unname(sqrt(phi2c / min(kc - 1, rc - 1))))
  expect_match(adj$effect_size$name, "bias-corrected")
  expect_lt(adj$effect_size$estimate, raw$effect_size$estimate)
  expect_equal(raw$effect_size$name, "Cramer's V")
})
