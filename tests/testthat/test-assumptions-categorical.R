test_that("categorical workflows choose chi-square or Fisher based on expected counts", {
  big <- tibble::tibble(x = rep(c("a", "b"), each = 30), y = rep(c("yes", "no"), 30))
  small <- tibble::tibble(x = c("a", "a", "a", "b"), y = c("yes", "yes", "no", "no"))

  x1 <- test_categorical(big, x, y, plot = FALSE)
  expect_equal(x1$recommended$test, "Chi-square test of independence")
  expect_true(any(grepl("Expected cell counts", x1$assumptions$name)))

  x2 <- test_categorical(small, x, y, plot = FALSE)
  expect_equal(x2$recommended$test, "Fisher exact test")
})

test_that("the expected-counts assumption panel agrees with the chi-square/Fisher recommendation", {
  # 10% of expected cells below 5 (satisfies the old Cochran 80% rule) but
  # still triggers Fisher under the any-cell-below-threshold rule actually
  # used for the recommendation; the assumption panel must not contradict it.
  tab <- matrix(c(200, 190, 180, 170, 8, 40, 38, 36, 34, 2), nrow = 2, byrow = TRUE)
  dimnames(tab) <- list(row = c("A", "B"), col = paste0("c", 1:5))
  long <- as.data.frame(as.table(tab))
  dat <- long[rep(seq_len(nrow(long)), long$Freq), c("row", "col")]

  x <- test_categorical(row ~ col, data = dat, plot = FALSE)
  expect_equal(x$recommended$test, "Fisher exact test")
  expected_panel <- x$assumptions[x$assumptions$name == "Expected cell counts", ]
  expect_equal(expected_panel$status, "not acceptable")
})

test_that("paired categorical and repeated categorical report assumption checks", {
  dat <- make_assumption_data(60)
  x <- test_paired_categorical(dat, controlled_baseline, controlled_3m, plot = FALSE)
  y <- test_repeated_categorical(dat, c(controlled_baseline, controlled_3m, controlled_6m), plot = FALSE)
  expect_true(nrow(format_assumptions(x$assumptions)) >= 1)
  expect_true(nrow(format_assumptions(y$assumptions)) >= 1)
})
