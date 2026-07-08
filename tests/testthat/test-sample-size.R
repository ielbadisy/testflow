test_that("continuous sample size supports paired and parallel planning", {
  x <- sample_size_continuous(
    design = "parallel",
    objective = "superiority",
    delta = 5,
    sd = 10,
    alpha = 0.05,
    power = 0.90
  )

  expect_s3_class(x, "sample_size")
  expect_equal(unname(x$n_adjusted["A"]), 85)
  expect_equal(unname(x$n_adjusted["B"]), 85)
  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
  tbl <- as_tibble(x)
  expect_equal(nrow(tbl), 1)
  expect_true(all(c("endpoint", "design", "objective", "formula", "reference") %in% names(tbl)))

  y <- sample_size_continuous(
    design = "paired",
    objective = "superiority",
    delta = 5,
    sd_diff = 10,
    alpha = 0.05,
    power = 0.90
  )
  expect_equal(y$n_adjusted, 43)
  expect_true(grepl("paired normal approximation", y$method, fixed = TRUE))

  r <- sample_size_continuous(
    design = "repeated",
    n_time = 4,
    correlation = 0.5,
    objective = "superiority",
    delta = 5,
    sd_diff = 10,
    alpha = 0.05,
    power = 0.90
  )
  expect_equal(r$design, "repeated (4 time points)")
  expect_true(grepl("repeated-measures normal approximation", r$method, fixed = TRUE))
  expect_equal(r$n_adjusted, 27)
  expect_true(inherits(plot(r, type = "summary"), "ggplot"))
})

test_that("binary paired planning uses discordant pairs", {
  x <- sample_size_binary(
    design = "paired",
    objective = "superiority",
    p10 = 0.2,
    p01 = 0.1,
    alpha = 0.05,
    power = 0.90
  )

  expect_s3_class(x, "sample_size")
  expect_equal(x$n_adjusted, 316)
  expect_s3_class(plot(x), "ggplot")
  expect_true(grepl("discordant pairs", x$method, fixed = TRUE))
})

test_that("survival and ordinal helpers return planned counts", {
  s <- sample_size_survival(
    hr = 0.7,
    survival_a = 0.8,
    survival_b = 0.7,
    alpha = 0.05,
    power = 0.90
  )
  expect_s3_class(s, "sample_size")
  expect_true(s$required_events > 0)
  expect_true(s$n_adjusted > 0)

  o <- sample_size_ordinal(
    p_superiority = 0.6,
    alpha = 0.05,
    power = 0.90
  )
  expect_equal(unname(o$n_adjusted["A"]), 176)
  expect_equal(o$n_total, 352)
})

test_that("sample_size dispatches to endpoint-specific helpers", {
  x <- sample_size(
    endpoint = "continuous",
    design = "paired",
    objective = "superiority",
    delta = 5,
    sd_diff = 10,
    alpha = 0.05,
    power = 0.90
  )

  expect_s3_class(x, "sample_size")
  expect_equal(x$n_adjusted, 43)
})
