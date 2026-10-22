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

test_that("Wilson/exact search returns a stable n: target met at n and beyond, not at n - 1", {
  scenarios <- list(
    list(p = 0.30, width = 0.05),
    list(p = 0.05, width = 0.02),
    list(p = 0.01, width = 0.005),
    list(p = 0.001, width = 0.001)
  )
  for (s in scenarios) {
    for (method in c("wilson", "exact")) {
      n <- testflow:::sample_size_binomial_required_n(s$p, s$width, 0.05, method, "expected", 1e7)
      after <- n:(n + max(50, ceiling(4 / s$p)))
      hw_after <- vapply(after, testflow:::sample_size_binomial_max_half_width, numeric(1),
                         p = s$p, alpha = 0.05, method = method, criterion = "expected")
      expect_true(all(hw_after <= s$width), info = paste(method, s$p, s$width))
      if (n > 2) {
        hw_n1 <- testflow:::sample_size_binomial_max_half_width(n - 1, s$p, 0.05, method, "expected")
        expect_true(hw_n1 > s$width, info = paste(method, s$p, s$width))
      }
    }
  }
})

test_that("Wilson expected criterion is evaluated at p and matches the closed-form interval", {
  wilson_gap <- function(p, n, z = stats::qnorm(0.975)) {
    centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    marge <- z / (1 + z^2 / n) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
    max(p - (centre - marge), centre + marge - p)
  }
  p <- 27 / 5418
  n <- testflow:::sample_size_binomial_required_n(p, 0.0015, 0.05, "wilson", "expected", 1e7)
  expect_equal(n, 10998L)
  expect_true(wilson_gap(p, n) <= 0.0015)
  expect_true(wilson_gap(p, n - 1) > 0.0015)
  expect_equal(testflow:::sample_size_binomial_max_half_width(n, p, 0.05, "wilson", "expected"),
               wilson_gap(p, n), tolerance = 1e-12)
  expect_equal(testflow:::sample_size_binomial_required_n(55 / 5419, 0.0025, 0.05, "wilson", "expected", 1e7), 7677L)

  # monotone in n: no rounding noise
  hw <- vapply(10000:12000, testflow:::sample_size_binomial_max_half_width, numeric(1),
               p = p, alpha = 0.05, method = "wilson", criterion = "expected")
  expect_true(all(diff(hw) < 0))

  x <- suppressWarnings(sample_size_precision(endpoint = "binary", design = "one_sample", p = p,
                                              width = 0.0015, method = "wilson"))
  expect_equal(x$diagnostics$complete_case_n, 10998)
  expect_equal(x$diagnostics$achieved_maximum_half_width, wilson_gap(p, 10998), tolerance = 1e-12)
})

test_that("worst_case criterion is at least as large as expected criterion for the same n", {
  n <- 500
  hw_expected <- testflow:::sample_size_binomial_max_half_width(n, 0.05, 0.05, "wilson", "expected")
  hw_worst <- testflow:::sample_size_binomial_max_half_width(n, 0.05, 0.05, "wilson", "worst_case")
  expect_true(hw_worst >= hw_expected)
})

test_that("the search matches a brute-force scan for the stable n under non-monotonicity", {
  # With integer counts (exact, or worst_case), max_half_width(n) is not
  # monotone: round(n * p) and the qbinom() band jump, so the target can be
  # met at n and missed again at a larger n. The search must return the
  # smallest n after which no failure occurs. Widths sit exactly at known
  # local dips (p = 0.05, n = 2:400), where the first passing n is not stable.
  brute_force_stable_n <- function(p, width, alpha, method, criterion, n_max) {
    ok <- vapply(2:n_max, function(n) {
      testflow:::sample_size_binomial_max_half_width(n, p, alpha, method, criterion) <= width
    }, logical(1))
    last_fail <- max(c(1L, (2:n_max)[!ok]))
    last_fail + 1L
  }

  cases <- list(
    list(p = 0.05, method = "exact", criterion = "expected", n0 = 70),
    list(p = 0.05, method = "wilson", criterion = "worst_case", n0 = 43),
    list(p = 0.05, method = "exact", criterion = "worst_case", n0 = 33),
    list(p = 0.01, method = "wilson", criterion = "worst_case", n0 = 20),
    list(p = 27 / 5418, method = "exact", criterion = "expected", n0 = 11000)
  )
  for (cfg in cases) {
    hw0 <- testflow:::sample_size_binomial_max_half_width(cfg$n0, cfg$p, 0.05, cfg$method, cfg$criterion)
    width <- hw0 + 1e-6
    n_search <- testflow:::sample_size_binomial_required_n(cfg$p, width, 0.05, cfg$method, cfg$criterion, 1e5)
    n_brute <- brute_force_stable_n(cfg$p, width, 0.05, cfg$method, cfg$criterion,
                                    n_search + 3 * max(50, ceiling(4 / cfg$p)))
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
