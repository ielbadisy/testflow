test_that("print and as_tibble expose core text", {
  dat <- make_cardio_data(80)
  x <- test_two_groups(dat, sbp_3m, sex)
  txt <- capture.output(print(x))
  expect_true(any(grepl("Two Independent Groups", txt, fixed = TRUE)))
  expect_true(any(grepl("sbp_3m", txt)))
  expect_true(any(grepl(x$recommended$test, txt, fixed = TRUE)))
  expect_true(any(grepl("Report", txt)))
  expect_true(any(grepl("H0:", txt, fixed = TRUE)))
  tbl <- as_tibble(x)
  expect_equal(nrow(tbl), 1)
  expect_true("null_hypothesis" %in% names(tbl))
  expect_true("conf.low" %in% names(tbl))
})

test_that("summary.testflow prints a compact vertical result", {
  dat <- make_cardio_data(80)
  x <- test_two_groups(dat, sbp_3m, sex)

  txt <- capture.output(summary(x))

  expect_true(any(grepl("Two Independent Groups (summary)", txt, fixed = TRUE)))
  expect_true(any(grepl("Workflow:", txt, fixed = TRUE)))
  expect_true(any(grepl("Recommended test:", txt, fixed = TRUE)))
  expect_true(any(grepl(x$recommended$test, txt, fixed = TRUE)))
  expect_true(any(grepl("H0:", txt, fixed = TRUE)))
  expect_true(any(grepl("p:", txt, fixed = TRUE)))
  expect_true(any(grepl("Effect size:", txt, fixed = TRUE)))
  expect_true(any(grepl("Report", txt, fixed = TRUE)))
})

test_that("testflow cli colors can be forced on and disabled", {
  dat <- make_cardio_data(80)
  x <- test_two_groups(dat, sbp_3m, sex)

  old <- options(testflow.cli_colors = TRUE, cli.num_colors = 1)
  on.exit(options(old), add = TRUE)
  colored <- capture.output(print(x))
  expect_true(any(cli::ansi_has_any(colored)))

  options(testflow.cli_colors = FALSE, cli.num_colors = 256)
  plain <- capture.output(print(x))
  expect_false(any(cli::ansi_has_any(plain)))
})

test_that("assumption labels do not print as NA", {
  dat <- make_cardio_data(80)
  one <- test_one_sample(dat, sbp_3m, mu = 140)
  two <- test_two_groups(dat, sbp_3m, sex)

  one_txt <- capture.output(print(one))
  two_txt <- capture.output(print(two))

  expect_false(any(grepl("^\\* NA:", one_txt)))
  expect_false(any(grepl("^\\* NA:", two_txt)))
  expect_true(any(grepl("Normality:", one_txt, fixed = TRUE)))
  expect_true(any(grepl("Normality:", two_txt, fixed = TRUE)))
})

test_that("screening workflows do not print fake test statistics", {
  dat <- make_cardio_data(80)

  corr <- suppressWarnings(test_correlation_matrix(dat, c(age, sbp_3m, ldl), plot = FALSE))
  outliers <- suppressWarnings(test_outliers(c(sbp_3m, ldl), data = dat, plot = FALSE))

  corr_txt <- capture.output(print(corr))
  outlier_txt <- capture.output(print(outliers))

  expect_false(any(grepl("statistic = NA", corr_txt, fixed = TRUE)))
  expect_false(any(grepl("p = NA", outlier_txt, fixed = TRUE)))
  expect_true(any(grepl("pairwise correlations reported", corr_txt, fixed = TRUE)))
  expect_true(any(grepl("flagged rows", outlier_txt, fixed = TRUE)))
})

test_that("print.testflow shows a workflow-specific title and relabeled fields", {
  set.seed(1)
  n <- 60
  dat <- tibble::tibble(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
  dat$y <- 1 + 0.5 * dat$x1 + rnorm(n, sd = 0.5)
  reg <- test_linear_regression(y ~ x1 + x2, data = dat)
  reg_txt <- capture.output(print(reg))
  expect_true(any(grepl("^Linear Regression$", reg_txt)))
  expect_true(any(grepl("Predictors: x1, x2", reg_txt, fixed = TRUE)))
  expect_false(any(grepl("^Group:", reg_txt)))

  icc_dat <- tibble::tibble(r1 = rnorm(20, 50, 10), r2 = rnorm(20, 50, 10), r3 = rnorm(20, 50, 10))
  icc <- test_icc(icc_dat, c(r1, r2, r3))
  icc_txt <- capture.output(print(icc))
  expect_true(any(grepl("^Intraclass Correlation$", icc_txt)))
  expect_false(any(grepl("^Group:", icc_txt)))

  reg_summary_txt <- capture.output(print(summary(reg)))
  expect_true(any(grepl("Linear Regression (summary)", reg_summary_txt, fixed = TRUE)))
  expect_true(any(grepl("Predictors:", reg_summary_txt, fixed = TRUE)))
})

test_that("print.testflow surfaces per-term results tables", {
  set.seed(1)
  n <- 80
  dat <- tibble::tibble(y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n))
  dat$y <- 2 + 0.6 * dat$x1 - 0.3 * dat$x2 + rnorm(n, sd = 0.8)
  reg <- test_linear_regression(y ~ x1 + x2, data = dat)
  reg_txt <- capture.output(print(reg))
  expect_true(any(grepl("Coefficients", reg_txt, fixed = TRUE)))
  expect_true(any(grepl("x1", reg_txt, fixed = TRUE)) && any(grepl("x2", reg_txt, fixed = TRUE)))
  # p-values in the table use the same "<0.001" convention as the primary
  # result, not a truncated "0.00" from the generic 2-decimal formatter.
  expect_true(any(grepl("<0.001", reg_txt, fixed = TRUE)))

  ns <- 150
  sd <- tibble::tibble(time = rexp(ns, 0.08), status = rbinom(ns, 1, 0.75), age = rnorm(ns, 60, 10))
  cox <- test_cox(Surv(time, status) ~ age, data = sd)
  cox_txt <- capture.output(print(cox))
  expect_true(any(grepl("Hazard ratios", cox_txt, fixed = TRUE)))

  tp <- 45; fp <- 10; fn <- 8; tn <- 90
  diag_dat <- tibble::tibble(
    test = c(rep("positive", tp + fp), rep("negative", fn + tn)),
    reference = c(rep("positive", tp), rep("negative", fp), rep("positive", fn), rep("negative", tn))
  )
  diag <- test_diagnostic(diag_dat, test, reference)
  diag_txt <- capture.output(print(diag))
  expect_true(any(grepl("Diagnostic accuracy", diag_txt, fixed = TRUE)))
  expect_true(any(grepl("Sensitivity", diag_txt, fixed = TRUE)))

  # A wide table (long text column + several numeric columns) must not wrap
  # onto a second block of misaligned rows.
  icc_dat <- tibble::tibble(r1 = rnorm(20, 50, 10), r2 = rnorm(20, 50, 10), r3 = rnorm(20, 50, 10))
  icc_txt <- capture.output(print(test_icc(icc_dat, c(r1, r2, r3))))
  icc_header_line <- icc_txt[which(grepl("method", icc_txt, fixed = TRUE))]
  expect_true(any(grepl("p.value", icc_header_line, fixed = TRUE)))

  # A matrix/table-sourced result gets meaningful column names, not the
  # generic Var1/Var2/Freq produced by as.data.frame(as.table(...)).
  agree_dat <- tibble::tibble(rater1 = sample(c("A", "B"), 100, replace = TRUE), rater2 = sample(c("A", "B"), 100, replace = TRUE))
  agree_txt <- capture.output(print(test_agreement(agree_dat, rater1, rater2)))
  expect_true(any(grepl("Rater 1", agree_txt, fixed = TRUE)) && any(grepl("Rater 2", agree_txt, fixed = TRUE)))
  expect_false(any(grepl("Var1", agree_txt, fixed = TRUE)))

  # groups' post-hoc is a wrapped TukeyHSD/pairwise.htest object, not a
  # tidy data frame - must be silently skipped, not error or dump raw output.
  cardio <- make_cardio_data(90)
  grp_txt <- capture.output(print(test_groups(sbp_3m ~ treatment, data = cardio)))
  expect_false(any(grepl("Post-hoc comparisons", grp_txt, fixed = TRUE)))

  # repeated's post-hoc IS a tidy tibble and should render.
  rep_txt <- capture.output(print(test_repeated(cardio, c(sbp_baseline, sbp_3m, sbp_6m), id = id)))
  expect_true(any(grepl("Post-hoc comparisons", rep_txt, fixed = TRUE)))
})

test_that("wide repeated workflow prints measure labels", {
  dat <- make_cardio_data(80)
  x <- test_repeated(dat, c(sbp_baseline, sbp_3m, sbp_6m), id = id)

  txt <- capture.output(print(x))

  expect_true(any(grepl("Outcome: sbp_baseline, sbp_3m, sbp_6m", txt, fixed = TRUE)))
  expect_false(any(grepl("Normality: value", txt, fixed = TRUE)))
  expect_false(any(grepl("population mean or location of value", txt, fixed = TRUE)))
})
