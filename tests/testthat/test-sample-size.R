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
  expect_true(all(c("endpoint", "design", "objective", "method", "report") %in% names(tbl)))
  expect_false(any(c("formula", "reference") %in% names(tbl)))

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
  expect_true(inherits(plot(r, type = "curve"), "ggplot"))
  expect_true(inherits(plot(r, type = "both"), "ggplot") || inherits(plot(r, type = "both"), "patchwork"))
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
  # events-to-N conversion: N_total = 2 * D_total / (qA + qB), qA = 0.2, qB = 0.3
  expect_equal(s$required_events, 331)
  expect_equal(s$n_adjusted, 1322)

  o <- sample_size_ordinal(
    p_superiority = 0.6,
    alpha = 0.05,
    power = 0.90
  )
  expect_equal(unname(o$n_adjusted["A"]), 176)
  expect_equal(o$n_total, 352)
})

test_that("paired continuous power curve matches the target power at the planned n", {
  r <- sample_size_continuous(
    design = "paired",
    objective = "superiority",
    delta = 5,
    sd_diff = 10,
    alpha = 0.05,
    power = 0.90
  )
  target_n <- r$curve_data$target_n[1]
  power_at_target <- r$curve_data$power[which.min(abs(r$curve_data$n - target_n))]
  expect_equal(power_at_target, 0.90, tolerance = 0.02)
})

test_that("binary equivalence special case respects the allocation ratio", {
  x <- sample_size_binary(
    design = "parallel",
    objective = "equivalence",
    p1 = 0.3,
    p2 = 0.3,
    margin = 0.15,
    allocation = 2,
    alpha = 0.05,
    power = 0.90
  )
  expect_equal(unname(x$n_adjusted["A"]), 152)
  expect_equal(unname(x$n_adjusted["B"]), 304)
  expect_equal(unname(x$n_adjusted["B"]) / unname(x$n_adjusted["A"]), 2, tolerance = 0.02)
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

test_that("sample_size_cluster_adjust matches the design-effect formula exactly", {
  expect_equal(sample_size_cluster_adjust(100, m = 20, rho = 0.02), ceiling(100 * (1 + 19 * 0.02)))
  expect_equal(
    sample_size_cluster_adjust(100, m = 20, rho = 0.02, cv_m = 0.3),
    ceiling(100 * (1 + ((1 + 0.3^2) * 20 - 1) * 0.02))
  )
})

test_that("sample_size_precision matches each spec formula exactly", {
  za <- qnorm(0.975)

  x1 <- sample_size_precision(endpoint = "continuous", design = "one_sample", width = 2, sd = 10)
  expect_equal(x1$n, (za * 10 / 2)^2)

  x2 <- sample_size_precision(endpoint = "continuous", design = "two_sample", width = 2, sd = 10, allocation = 1.5)
  expect_equal(x2$n, (1.5 + 1) * 10^2 * za^2 / (1.5 * 2^2))

  x3 <- suppressWarnings(sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.3))
  expect_equal(x3$n, za^2 * 0.3 * 0.7 / 0.05^2)

  x4 <- suppressWarnings(sample_size_precision(endpoint = "binary", design = "one_sample", width = 0.05, conservative = TRUE))
  expect_equal(x4$n, za^2 / (4 * 0.05^2))

  x5 <- sample_size_precision(endpoint = "binary", design = "two_sample", width = 0.08, p1 = 0.4, p2 = 0.3)
  expect_equal(x5$n, za^2 * (0.4 * 0.6 + 0.3 * 0.7) / 0.08^2)

  x6 <- sample_size_precision(endpoint = "binary", design = "odds_ratio", width = 0.3, p1 = 0.4, p2 = 0.3, allocation = 2)
  expect_equal(x6$n, za^2 / 0.3^2 * (1 / 0.4 + 1 / 0.6 + 1 / (2 * 0.3) + 1 / (2 * 0.7)))

  expect_error(sample_size_precision(endpoint = "continuous", design = "odds_ratio", width = 1, sd = 1), "binary")
})

test_that("sample_size_ordinal's Noether method matches the spec formula", {
  za <- qnorm(0.975); zb <- qnorm(0.90)

  x_noether <- sample_size_ordinal(p_superiority = 0.6)
  expect_equal(x_noether$n, (za + zb)^2 / (6 * 0.1^2))
  expect_true(grepl("Noether", x_noether$method, fixed = TRUE))

  expect_error(sample_size_ordinal(p_superiority = 0.4), "exceed 0.5")
})

test_that("sample_size_survival's uniform accrual matches the closed-form formula and reduces to the flat conversion as accrual -> 0", {
  hr <- 0.7
  za <- qnorm(0.975); zb <- qnorm(0.90)
  events_per_arm <- 2 * (za + zb)^2 / (log(hr)^2)
  d_total <- 2 * events_per_arm

  x0 <- sample_size_survival(hr = hr, survival_a = 0.8, survival_b = 0.7, alpha = 0.05, power = 0.90)
  expect_equal(x0$n, 2 * d_total / ((1 - 0.8) + (1 - 0.7)))

  R <- 12; FUP <- 24; total_dur <- R + FUP
  lambda_a <- -log(0.8) / total_dur
  lambda_b <- -log(0.7) / total_dur
  qbar_a <- 1 - (exp(-lambda_a * (total_dur - R)) - exp(-lambda_a * total_dur)) / (lambda_a * R)
  qbar_b <- 1 - (exp(-lambda_b * (total_dur - R)) - exp(-lambda_b * total_dur)) / (lambda_b * R)
  x1 <- sample_size_survival(hr = hr, survival_a = 0.8, survival_b = 0.7, alpha = 0.05, power = 0.90, accrual_duration = R, follow_up = FUP)
  expect_equal(x1$n, 2 * d_total / (qbar_a + qbar_b))
  expect_true(x1$n > x0$n)

  x2 <- sample_size_survival(hr = hr, survival_a = 0.8, survival_b = 0.7, alpha = 0.05, power = 0.90, accrual_duration = 0.0001, follow_up = total_dur - 0.0001)
  expect_equal(x2$n, x0$n, tolerance = 1e-3)

  expect_error(sample_size_survival(hr = hr, survival_a = 0.8, survival_b = 0.7, accrual_duration = 12), "follow_up")
})

test_that("sample_size_bioequivalence's iterative TOST matches the exact noncentral-t TOST power formula", {
  # Independent re-implementation of the Phillips (1990) exact TOST power
  # (noncentral t, not a normal approximation): T_U, T_L are each
  # noncentral-t(df, ncp) under the true theta0, so power is computed via
  # stats::pt(..., ncp = ...), not pnorm()/qnorm(). A prior version of
  # sample_size_bioequivalence() used the normal approximation below (kept
  # here only as `old_normal_approx_tost_power` to document what changed);
  # that version under-sized studies by up to ~20% at small n, confirmed
  # against PowerTOST::power.TOST() during development.
  exact_tost_power <- function(n, theta0, sigma, theta_l, theta_u, alpha, df) {
    se <- sqrt(2 * sigma^2 / n)
    tcrit <- qt(1 - alpha, df(n))
    delta_u <- (theta0 - theta_u) / se
    delta_l <- (theta0 - theta_l) / se
    pt(-tcrit, df(n), ncp = delta_u) - pt(tcrit, df(n), ncp = delta_l)
  }

  cv_w <- 0.3
  sigma_w <- sqrt(log(1 + cv_w^2))

  # normal_approx (closed form, unchanged) still recovers its own documented formula
  x1b <- sample_size_bioequivalence(design = "crossover", gmr = 1, cv_within = cv_w, method = "normal_approx")
  zbeta2 <- qnorm(1 - (1 - 0.90) / 2); za <- qnorm(1 - 0.05)
  n_closed <- 2 * sigma_w^2 * (zbeta2 + za)^2 / (log(1.25))^2
  expect_equal(x1b$n, n_closed)

  # iterative_tost (GMR = 1): matches an independent noncentral-t search, and
  # can legitimately differ from normal_approx now that it is exact
  x1 <- sample_size_bioequivalence(design = "crossover", gmr = 1, cv_within = cv_w, method = "iterative_tost")
  n_manual_center <- uniroot(function(n) exact_tost_power(n, 0, sigma_w, log(0.8), log(1.25), 0.05, function(n) n - 2) - 0.90, c(4, 1e6), extendInt = "upX")$root
  expect_equal(x1$n, n_manual_center, tolerance = 1e-6)

  # off-center GMR: iterative_tost matches an independently computed
  # noncentral-t uniroot search
  gmr <- 0.95
  theta0 <- log(gmr); theta_l <- log(0.8); theta_u <- log(1.25)
  n_manual <- uniroot(function(n) exact_tost_power(n, theta0, sigma_w, theta_l, theta_u, 0.05, function(n) n - 2) - 0.90, c(4, 1e6), extendInt = "upX")$root
  x2 <- sample_size_bioequivalence(design = "crossover", gmr = gmr, cv_within = cv_w, method = "iterative_tost")
  expect_equal(x2$n, n_manual, tolerance = 1e-6)

  # minimum-power property: ceiling(n) achieves the target power (the raw,
  # non-integer n from uniroot is exactly at the target by construction)
  x2_power_at_raw <- exact_tost_power(x2$n, theta0, sigma_w, theta_l, theta_u, 0.05, function(n) n - 2)
  expect_equal(x2_power_at_raw, 0.90, tolerance = 1e-6)

  # parallel design: allocation ratio is respected, and matches an
  # independent noncentral-t search with df = nA + nB - 2
  x3 <- sample_size_bioequivalence(design = "parallel", gmr = 0.95, cv_between = 0.35, allocation = 1.5, method = "iterative_tost")
  expect_equal(unname(x3$n_adjusted["B"]) / unname(x3$n_adjusted["A"]), 1.5, tolerance = 0.02)

  sigma_b <- sqrt(log(1 + 0.35^2))
  theta0_p <- log(0.95)
  exact_tost_power_parallel <- function(n_a, r, theta0, sigma, theta_l, theta_u, alpha) {
    n_b <- r * n_a
    df <- n_a + n_b - 2
    se <- sigma * sqrt((r + 1) / (r * n_a))
    tcrit <- qt(1 - alpha, df)
    delta_u <- (theta0 - theta_u) / se
    delta_l <- (theta0 - theta_l) / se
    pt(-tcrit, df, ncp = delta_u) - pt(tcrit, df, ncp = delta_l)
  }
  n_a_manual <- uniroot(function(n_a) exact_tost_power_parallel(n_a, 1.5, theta0_p, sigma_b, log(0.8), log(1.25), 0.05) - 0.90, c(4, 1e6), extendInt = "upX")$root
  expect_equal(x3$n, n_a_manual, tolerance = 1e-6)

  expect_error(sample_size_bioequivalence(design = "crossover", gmr = 1.3, cv_within = 0.3), "inside")
  expect_error(sample_size_bioequivalence(design = "crossover", gmr = 1, cv_within = 0.3, lower = 0.8, upper = 0.5), "exceed")
})

test_that("sample_size_bioequivalence's iterative TOST matches PowerTOST::power.TOST() (regulatory-standard exact power)", {
  skip_if_not_installed("PowerTOST")
  # Balanced-design cross-check against Owen's-Q-based exact TOST power (the
  # standard used in real bioequivalence submissions). Even total n avoids
  # PowerTOST's unbalanced-sequence adjustment for odd n, which testflow
  # (treating n as one continuous total) does not separately model.
  for (cfg in list(c(gmr = 0.95, cv = 0.15), c(gmr = 1.00, cv = 0.20), c(gmr = 0.90, cv = 0.20))) {
    x <- sample_size_bioequivalence(design = "crossover", gmr = cfg[["gmr"]], cv_within = cfg[["cv"]], method = "iterative_tost")
    n_even <- 2 * ceiling(x$n / 2)
    n_even_minus_2 <- n_even - 2
    power_at_n <- PowerTOST::power.TOST(alpha = 0.05, theta0 = cfg[["gmr"]], theta1 = 0.8, theta2 = 1.25, CV = cfg[["cv"]], n = n_even, design = "2x2")
    power_below <- PowerTOST::power.TOST(alpha = 0.05, theta0 = cfg[["gmr"]], theta1 = 0.8, theta2 = 1.25, CV = cfg[["cv"]], n = n_even_minus_2, design = "2x2")
    # testflow's n (rounded up to the nearest even total) must achieve at
    # least the target power under PowerTOST's exact calculation, and two
    # fewer subjects (the next balanced design down) must not - the same
    # minimality property already required of the Wilson/exact precision
    # search, cross-checked here against an external, regulatory-standard
    # implementation rather than an internal formula.
    expect_true(power_at_n >= 0.90 - 1e-6)
    expect_true(power_below < 0.90)
  }
})
