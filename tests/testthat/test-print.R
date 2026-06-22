test_that("print and as_tibble expose core text", {
  dat <- make_cardio_data(80)
  x <- test_two_groups(dat, sbp_3m, sex)
  txt <- capture.output(print(x))
  expect_true(any(grepl("Statistical test workflow", txt)))
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

  expect_true(any(grepl("testflow summary", txt, fixed = TRUE)))
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
