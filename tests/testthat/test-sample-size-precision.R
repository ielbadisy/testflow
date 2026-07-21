za <- qnorm(0.975)

# ---- Backward compatibility -------------------------------------------------

test_that("Wald one-proportion precision matches the closed-form formula exactly", {
  x <- suppressWarnings(sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.30))
  expect_equal(x$n, za^2 * 0.3 * 0.7 / 0.05^2)
  expect_equal(x$diagnostics$ci_method, "wald")
  expect_equal(x$diagnostics$precision_criterion, "expected")
})

test_that("Conservative Wald one-proportion precision equals z^2 / (4w^2)", {
  x <- suppressWarnings(sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, conservative = TRUE))
  expect_equal(x$n, za^2 / (4 * 0.05^2))
  expect_equal(x$diagnostics$anticipated_prevalence, 0.5)
})

test_that("Existing continuous precision results remain unchanged", {
  x1 <- sample_size_precision(endpoint = "continuous", design = "one_sample", width = 2, sd = 10)
  expect_equal(x1$n, (za * 10 / 2)^2)
  x2 <- sample_size_precision(endpoint = "continuous", design = "two_sample", width = 2, sd = 10, allocation = 1.5)
  expect_equal(x2$n, (1.5 + 1) * 10^2 * za^2 / (1.5 * 2^2))
})

test_that("Existing odds-ratio precision results remain unchanged", {
  x <- sample_size_precision(endpoint = "binary", design = "odds_ratio", width = 0.3, p1 = 0.4, p2 = 0.3, allocation = 2)
  expect_equal(x$n, za^2 / 0.3^2 * (1 / 0.4 + 1 / 0.6 + 1 / (2 * 0.3) + 1 / (2 * 0.7)))
})

# ---- Wilson interval helper --------------------------------------------------

test_that("sample_size_binomial_ci(wilson) matches the independent Wilson formula", {
  z <- qnorm(0.975)
  cases <- list(c(0, 100), c(1, 100), c(30, 100), c(100, 100))
  for (cfg in cases) {
    x <- cfg[1]; n <- cfg[2]
    p_hat <- x / n
    denom <- 1 + z^2 / n
    center <- (p_hat + z^2 / (2 * n)) / denom
    half <- (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
    expected_lower <- max(0, min(1, center - half))
    expected_upper <- max(0, min(1, center + half))
    ci <- testflow:::sample_size_binomial_ci(x, n, 0.05, "wilson")
    expect_equal(ci$lower, expected_lower, tolerance = 1e-10)
    expect_equal(ci$upper, expected_upper, tolerance = 1e-10)
  }
})

# ---- Exact interval helper ---------------------------------------------------

test_that("sample_size_binomial_ci(exact) matches stats::binom.test()", {
  cases <- list(c(0, 50), c(1, 100), c(30, 100), c(50, 50))
  for (cfg in cases) {
    x <- cfg[1]; n <- cfg[2]
    bt <- stats::binom.test(x, n, conf.level = 0.95)$conf.int
    ci <- testflow:::sample_size_binomial_ci(x, n, 0.05, "exact")
    expect_equal(ci$lower, bt[1], tolerance = 1e-6)
    expect_equal(ci$upper, bt[2], tolerance = 1e-6)
  }
})

# ---- Minimum-n property -------------------------------------------------------

test_that("Wilson/exact minimum-n search satisfies the precision target and is minimal", {
  scenarios <- list(
    list(p = 0.30, width = 0.05),
    list(p = 0.05, width = 0.02),
    list(p = 0.01, width = 0.005),
    list(p = 0.001, width = 0.001)
  )
  for (s in scenarios) {
    for (method in c("wilson", "exact")) {
      n <- testflow:::sample_size_binomial_required_n(s$p, s$width, 0.05, method, "expected", 1e7)
      hw_n <- testflow:::sample_size_binomial_max_half_width(n, s$p, 0.05, method, "expected")
      expect_true(hw_n <= s$width, info = paste(method, s$p, s$width))
      if (n > 2) {
        hw_n1 <- testflow:::sample_size_binomial_max_half_width(n - 1, s$p, 0.05, method, "expected")
        expect_true(hw_n1 > s$width, info = paste(method, s$p, s$width))
      }
    }
  }
})

test_that("worst_case criterion is at least as large as expected criterion for the same n", {
  n <- 500
  hw_expected <- testflow:::sample_size_binomial_max_half_width(n, 0.05, 0.05, "wilson", "expected")
  hw_worst <- testflow:::sample_size_binomial_max_half_width(n, 0.05, 0.05, "wilson", "worst_case")
  expect_true(hw_worst >= hw_expected)
})

test_that("the search is robust to local non-monotonicity in max_half_width(n)", {
  # max_half_width(n) is not strictly monotone in n: round(n * p) and (for
  # criterion = "worst_case") the qbinom() quantile band jump discretely, so
  # an isolated smaller n can pass even when nearby larger n fail. A search
  # that only checks whether n - 1 also passes can converge to a boundary
  # more than one step above the true minimum. Confirm the search matches an
  # exhaustive brute-force scan at width values chosen to sit exactly at
  # known local dips (found by scanning max_half_width(n) for n = 2:400 at
  # p = 0.05, where up to ~7% of n -> n+1 steps are non-monotone).
  brute_force_min_n <- function(p, width, alpha, method, criterion, n_max) {
    for (n in 2:n_max) {
      if (testflow:::sample_size_binomial_max_half_width(n, p, alpha, method, criterion) <= width) return(n)
    }
    NA_integer_
  }

  cases <- list(
    list(p = 0.05, method = "wilson", criterion = "expected", n0 = 70),
    list(p = 0.05, method = "wilson", criterion = "worst_case", n0 = 43),
    list(p = 0.05, method = "exact", criterion = "worst_case", n0 = 33),
    list(p = 0.01, method = "wilson", criterion = "worst_case", n0 = 20)
  )
  for (cfg in cases) {
    hw0 <- testflow:::sample_size_binomial_max_half_width(cfg$n0, cfg$p, 0.05, cfg$method, cfg$criterion)
    width <- hw0 + 1e-6
    n_search <- testflow:::sample_size_binomial_required_n(cfg$p, width, 0.05, cfg$method, cfg$criterion, 1e5)
    n_brute <- brute_force_min_n(cfg$p, width, 0.05, cfg$method, cfg$criterion, n_search + 100)
    expect_equal(n_search, n_brute, info = paste(cfg$method, cfg$criterion, cfg$p, width))
  }
})

# ---- Rare-event diagnostics ---------------------------------------------------

test_that("rare-event diagnostics are calculated correctly", {
  x <- sample_size_precision(
    endpoint = "binary", design = "one_sample", width = 0.02, p = 0.01, method = "wilson"
  )
  n <- x$diagnostics$adjusted_n
  expect_equal(x$diagnostics$expected_events, n * 0.01)
  expect_equal(x$diagnostics$probability_zero_events, (1 - 0.01)^n, tolerance = 1e-8)
  expect_equal(x$diagnostics$probability_at_least_one_event, 1 - (1 - 0.01)^n, tolerance = 1e-8)
  expect_true(!is.null(x$diagnostics$approximation_warning))
})

test_that("no rare-event warning fires when both expected cells exceed the threshold", {
  x <- suppressWarnings(sample_size_precision(
    endpoint = "binary", design = "one_sample", width = 0.05, p = 0.30, method = "wilson"
  ))
  expect_null(x$diagnostics$approximation_warning)
  expect_warning(
    sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.30, method = "wilson"),
    NA
  )
})

test_that("Wald method always attaches an approximation warning", {
  expect_warning(
    sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.30),
    "Wald"
  )
})

# ---- Dropout -------------------------------------------------------------------

test_that("complete_case_n reflects the statistical criterion and adjusted_n applies dropout", {
  x <- suppressWarnings(sample_size_precision(
    endpoint = "binary", design = "one_sample", width = 0.02, p = 0.05, method = "wilson", dropout = 0.10
  ))
  n0 <- testflow:::sample_size_binomial_required_n(0.05, 0.02, 0.05, "wilson", "expected", 1e7)
  expect_equal(x$diagnostics$complete_case_n, n0)
  expect_equal(x$diagnostics$adjusted_n, ceiling(n0 / (1 - 0.10)))
  expect_true(x$diagnostics$achieved_maximum_half_width <= 0.02)
})

# ---- Validation ------------------------------------------------------------

test_that("binary one-sample precision validation errors are informative", {
  expect_error(sample_size_precision(endpoint = "binary", design = "one_sample", width = -0.01, p = 0.3), "positive")
  expect_error(sample_size_precision(endpoint = "binary", design = "one_sample", width = 1, p = 0.3), "\\(0, 1\\)")
  expect_error(sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0), "\\(0, 1\\)")
  expect_error(sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 1), "\\(0, 1\\)")
  expect_error(
    sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.3, conservative = TRUE, method = "wilson"),
    "conservative"
  )
  expect_error(
    sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.3, conservative = TRUE, method = "exact"),
    "conservative"
  )
  expect_error(
    sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.3, method = "wilson", max_n = 1),
    "max_n"
  )
  expect_error(
    sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.3, method = "wilson", max_n = Inf),
    "max_n"
  )
  expect_error(
    sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.001, p = 0.001, method = "wilson", max_n = 50),
    "No sample size"
  )
  expect_error(
    sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.3, min_expected_events = -1),
    "non-negative"
  )
})

# ---- Unequal allocation (two-sample binary precision) -------------------------

test_that("two-sample binary precision implements unequal allocation correctly", {
  width <- 0.08; p1 <- 0.4; p2 <- 0.3
  for (r in c(1, 1.5, 2)) {
    x <- sample_size_precision(
      endpoint = "binary", design = "two_sample", width = width, p1 = p1, p2 = p2, allocation = r
    )
    n_a_expected <- za^2 / width^2 * (p1 * (1 - p1) + p2 * (1 - p2) / r)
    n_b_expected <- r * n_a_expected
    expect_equal(x$n, n_a_expected)
    expect_equal(unname(x$n_adjusted["A"]), sample_size_round(n_a_expected))
    expect_equal(unname(x$n_adjusted["B"]), sample_size_round(n_b_expected))
    expect_equal(x$n_total, unname(x$n_adjusted["A"]) + unname(x$n_adjusted["B"]))
    if (r > 1) {
      expect_equal(unname(x$n_adjusted["B"]) / unname(x$n_adjusted["A"]), r, tolerance = 0.05)
    }
  }
})

test_that("two-sample binary allocation = 1 reproduces the original equal-allocation formula", {
  x <- sample_size_precision(endpoint = "binary", design = "two_sample", width = 0.08, p1 = 0.4, p2 = 0.3)
  expect_equal(x$n, za^2 * (0.4 * 0.6 + 0.3 * 0.7) / 0.08^2)
})
