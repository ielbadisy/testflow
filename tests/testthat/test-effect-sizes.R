test_that("Cohen's d formulas match their implementations", {
  one_sample <- cohens_d_one_sample(c(1, 2, 3), mu = 1)
  expect_equal(one_sample$estimate, 1)

  dat <- .bt(
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

test_that("rank-biserial correlation shares the sign convention of Cohen's d", {
  dat <- .bt(
    y = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    g = rep(c("a", "b"), each = 6)
  )
  d <- cohens_d_independent(dat, "y", "g")$estimate
  r <- rank_biserial_two_groups(dat, "y", "g")$estimate
  expect_lt(d, 0)
  expect_lt(r, 0)
  expect_equal(r, -1)
  skip_if_not_installed("effectsize")
  set.seed(1)
  dat2 <- .bt(y = c(stats::rexp(30, 1), stats::rexp(30, 0.3)), g = rep(c("a", "b"), each = 30))
  ref <- effectsize::rank_biserial(dat2$y[dat2$g == "a"], dat2$y[dat2$g == "b"])$r_rank_biserial
  expect_equal(rank_biserial_two_groups(dat2, "y", "g")$estimate, ref)
})
