test_that("samplesizecontinuous(method = 't') matches stats::power.t.test", {
  r <- samplesizecontinuous(delta = 0.5, sd = 1, power = 0.8, method = "t")
  ref <- stats::power.t.test(delta = 0.5, sd = 1, power = 0.8)$n
  expect_equal(unname(r$n_per_group[["A"]]), as.integer(ceiling(ref)))
  expect_match(r$method, "noncentral t")

  p <- samplesizecontinuous(design = "paired", delta = 0.4, sd_diff = 1, power = 0.9, method = "t")
  ref_p <- stats::power.t.test(delta = 0.4, sd = 1, power = 0.9, type = "paired")$n
  expect_equal(p$n_adjusted, as.integer(ceiling(ref_p)))

  ni <- samplesizecontinuous(objective = "noninferiority", delta = 0.5, sd = 1, power = 0.9, alpha = 0.025, method = "t")
  ref_ni <- stats::power.t.test(delta = 0.5, sd = 1, power = 0.9, sig.level = 0.025, alternative = "one.sided")$n
  expect_equal(unname(ni$n_per_group[["A"]]), as.integer(ceiling(ref_ni)))
})

test_that("method = 't' equivalence agrees with PowerTOST", {
  skip_if_not_installed("PowerTOST")
  eq <- samplesizecontinuous(objective = "equivalence", delta = 1, sd = 1.5, expected_difference = 0.2, power = 0.8, method = "t")
  ref <- PowerTOST::sampleN.TOST(alpha = 0.05, targetpower = 0.8, logscale = FALSE, theta0 = 0.2,
                                 theta1 = -1, theta2 = 1, CV = 1.5, design = "parallel", print = FALSE)
  expect_equal(2 * unname(eq$n_per_group[["A"]]), ref[["Sample size"]])
})

test_that("default method stays the normal approximation", {
  r <- samplesizecontinuous(delta = 0.5, sd = 1, power = 0.8)
  expect_equal(unname(r$n_per_group[["A"]]), 63L)
  expect_false(grepl("noncentral", r$method))
})
