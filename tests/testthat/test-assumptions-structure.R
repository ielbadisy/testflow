test_that("workflows return non-empty assumptions and print Assumptions", {
  dat <- make_assumption_data(80)
  dat$grp3 <- factor(ceiling(3 * rank(dat$age, ties.method = "first") / length(dat$age)))
  objs <- list(
    testonesample(dat, sbp_3m, mu = 140, plot = FALSE),
    testtwogroups(dat, sbp_3m, sex, plot = FALSE),
    testpaired(sbp_3m ~ sbp_baseline, data = dat, plot = FALSE),
    testgroups(dat, sbp_3m, grp3, plot = FALSE),
    testfactorial(dat, sbp_3m ~ sex * treatment, plot = FALSE),
    testrepeated(dat, c(sbp_baseline, sbp_3m, sbp_6m), id = id, plot = FALSE),
    testcategorical(dat, treatment, controlled_3m, plot = FALSE),
    testpairedcategorical(dat, controlled_baseline, controlled_3m, plot = FALSE),
    testrepeatedcategorical(dat, c(controlled_baseline, controlled_3m, controlled_6m), plot = FALSE),
    testproportion(dat, controlled_3m, success = "yes", p = 0.5, plot = FALSE),
    suppressWarnings(testmultinomial(dat, treatment, plot = FALSE)),
    testcorrelation(sbp_3m ~ age, data = dat, plot = FALSE),
    suppressWarnings(testcorrelationmatrix(dat, c(age, sbp_3m, ldl), plot = FALSE)),
    suppressWarnings(testoutliers(c(sbp_3m, ldl, crp), data = dat, plot = FALSE))
  )

  for (x in objs) {
    expect_false(is.null(x$assumptions))
    expect_gt(nrow(format_assumptions(x$assumptions)), 0)
    txt <- capture.output(print(x))
    expect_true(any(grepl("Assumptions", txt, fixed = TRUE)))
    expect_false(any(grepl("No assumptions reported.", txt, fixed = TRUE)))
  }
})
