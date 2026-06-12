test_that("test_categorical builds table and computes Cramer's V", {
  dat <- tibble::tibble(x = rep(c("a", "b"), each = 30), y = rep(c("yes", "no"), 30))
  z <- test_categorical(dat, x, y)
  expect_s3_class(z, "testflow_categorical")
  expect_equal(z$recommended$test, "Chi-square test of independence")
  expect_equal(z$effect_size$name, "Cramer's V")
})

test_that("test_categorical uses Fisher with small expected counts", {
  dat <- tibble::tibble(x = c("a", "a", "a", "b"), y = c("yes", "yes", "no", "no"))
  z <- test_categorical(dat, x, y)
  expect_equal(z$recommended$test, "Fisher exact test")
})

test_that("test_multinomial converts table counts to numeric internally", {
  dat <- tibble::tibble(x = c("a", "a", "b", "c", "c", "c"))
  z <- test_multinomial(dat, x, p = c(1 / 3, 1 / 3, 1 / 3))
  expect_s3_class(z, "testflow_multinomial")
  expect_equal(z$posthoc$expected, rep(2, 3))
})
