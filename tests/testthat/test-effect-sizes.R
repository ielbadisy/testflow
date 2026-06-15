test_that("Cohen's d formulas match their implementations", {
  one_sample <- cohens_d_one_sample(c(1, 2, 3), mu = 1)
  expect_equal(one_sample$estimate, 1)

  dat <- tibble::tibble(
    y = c(1, 2, 3, 2, 4, 6),
    g = rep(c("a", "b"), each = 3)
  )
  independent <- cohens_d_independent(dat, "y", "g")
  pooled <- sqrt(((3 - 1) * stats::var(c(1, 2, 3)) + (3 - 1) * stats::var(c(2, 4, 6))) / (3 + 3 - 2))
  expect_equal(independent$estimate, (mean(c(1, 2, 3)) - mean(c(2, 4, 6))) / pooled)

  paired <- cohens_d_paired(before = c(1, 2, 3), after = c(2, 4, 6))
  diff <- c(2, 4, 6) - c(1, 2, 3)
  expect_equal(paired$estimate, mean(diff) / stats::sd(diff))
})
