test_that("categorical workflows choose chi-square or Fisher based on expected counts", {
  big <- tibble::tibble(x = rep(c("a", "b"), each = 30), y = rep(c("yes", "no"), 30))
  small <- tibble::tibble(x = c("a", "a", "a", "b"), y = c("yes", "yes", "no", "no"))

  x1 <- test_categorical(big, x, y, plot = FALSE)
  expect_equal(x1$recommended$test, "Chi-square test of independence")
  expect_true(any(grepl("Expected cell counts", x1$assumptions$name)))

  x2 <- test_categorical(small, x, y, plot = FALSE)
  expect_equal(x2$recommended$test, "Fisher exact test")
})

test_that("paired categorical and repeated categorical report assumption checks", {
  dat <- make_assumption_data(60)
  x <- test_paired_categorical(dat, controlled_baseline, controlled_3m, plot = FALSE)
  y <- test_repeated_categorical(dat, c(controlled_baseline, controlled_3m, controlled_6m), plot = FALSE)
  expect_true(nrow(format_assumptions(x$assumptions)) >= 1)
  expect_true(nrow(format_assumptions(y$assumptions)) >= 1)
})
