test_that("workflows return non-empty assumptions and print Assumptions", {
  dat <- make_assumption_data(80)
  dat$grp3 <- factor(dplyr::ntile(dat$age, 3))
  objs <- list(
    test_one_sample(dat, sbp_3m, mu = 140, plot = FALSE),
    test_two_groups(dat, sbp_3m, sex, plot = FALSE),
    test_paired(sbp_3m ~ sbp_baseline, data = dat, plot = FALSE),
    test_groups(dat, sbp_3m, grp3, plot = FALSE),
    test_factorial(dat, sbp_3m ~ sex * treatment, plot = FALSE),
    test_repeated(dat, c(sbp_baseline, sbp_3m, sbp_6m), id = id, plot = FALSE),
    test_categorical(dat, treatment, controlled_3m, plot = FALSE),
    test_paired_categorical(dat, controlled_baseline, controlled_3m, plot = FALSE),
    test_repeated_categorical(dat, c(controlled_baseline, controlled_3m, controlled_6m), plot = FALSE),
    test_proportion(dat, controlled_3m, success = "yes", p = 0.5, plot = FALSE),
    test_multinomial(dat, treatment, plot = FALSE),
    test_correlation(sbp_3m ~ age, data = dat, plot = FALSE),
    test_correlation_matrix(dat, c(age, sbp_3m, ldl), plot = FALSE),
    test_outliers(c(sbp_3m, ldl, crp), data = dat, plot = FALSE)
  )

  for (x in objs) {
    expect_false(is.null(x$assumptions))
    expect_gt(nrow(format_assumptions(x$assumptions)), 0)
    txt <- capture.output(print(x))
    expect_true(any(grepl("Assumptions", txt, fixed = TRUE)))
    expect_false(any(grepl("No assumptions reported.", txt, fixed = TRUE)))
  }
})
