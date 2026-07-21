# Developer-facing numerical validation for one-proportion precision planning
# (sample_size_precision(endpoint = "binary", design = "one_sample")).
#
# Not part of the package's automated test suite. Run interactively with:
#   Rscript inst/validation/validate-sample-size-precision.R
# (from a source checkout, with the package loaded via devtools::load_all()).
#
# For each (p, width, alpha, method) combination, this script:
#   1. calls testflow::sample_size_precision() for the planned n;
#   2. independently recomputes the achieved confidence interval at the
#      anticipated event count, using formulas written out here (not by
#      calling the package's internal helpers) - an independent check, not a
#      round-trip through the same code;
#   3. checks the minimum-n property (n passes, n - 1 fails);
#   4. reports expected events and the probability of zero events.
#
# The script stops (via stopifnot()) the moment any invariant fails.

if (requireNamespace("devtools", quietly = TRUE) && !"testflow" %in% loadedNamespaces()) {
  devtools::load_all(quiet = TRUE)
} else if (!requireNamespace("testflow", quietly = TRUE)) {
  stop("Load the testflow package (devtools::load_all()) before running this script.")
}

independent_ci <- function(x, n, alpha, method) {
  z <- stats::qnorm(1 - alpha / 2)
  p_hat <- x / n
  if (method == "wald") {
    se <- sqrt(p_hat * (1 - p_hat) / n)
    lower <- p_hat - z * se
    upper <- p_hat + z * se
  } else if (method == "wilson") {
    denom <- 1 + z^2 / n
    center <- (p_hat + z^2 / (2 * n)) / denom
    half <- (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
    lower <- center - half
    upper <- center + half
  } else if (method == "exact") {
    lower <- if (x == 0) 0 else stats::qbeta(alpha / 2, x, n - x + 1)
    upper <- if (x == n) 1 else stats::qbeta(1 - alpha / 2, x + 1, n - x)
  } else {
    stop("Unknown method: ", method)
  }
  c(lower = max(0, min(1, lower)), upper = max(0, min(1, upper)))
}

# For "wald", the reference proportion is p itself (the closed-form n was
# derived at continuous p); for "wilson"/"exact", x must be an integer count,
# so the reference is round(n * p) / n.
independent_max_half_width <- function(n, p, alpha, method) {
  if (method == "wald") {
    x <- p * n
    p_ref <- p
  } else {
    x <- round(n * p)
    p_ref <- x / n
  }
  ci <- independent_ci(x, n, alpha, method)
  max(p_ref - ci["lower"], ci["upper"] - p_ref)
}

grid <- expand.grid(
  p = c(0.001, 0.005, 0.01, 0.05, 0.10, 0.30, 0.50),
  alpha = c(0.01, 0.05, 0.10),
  method = c("wald", "wilson", "exact"),
  stringsAsFactors = FALSE
)
# Width scaled to each prevalence so the target is always plausible
# (roughly half the anticipated prevalence, floored for very small p).
grid$width <- pmax(0.15 * pmin(grid$p, 1 - grid$p), 0.002)

results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  p <- grid$p[i]; width <- grid$width[i]; alpha <- grid$alpha[i]; method <- grid$method[i]

  res <- suppressWarnings(testflow::sample_size_precision(
    endpoint = "binary", design = "one_sample",
    width = width, p = p, alpha = alpha, method = method, criterion = "expected"
  ))
  n <- res$diagnostics$complete_case_n

  achieved <- independent_max_half_width(n, p, alpha, method)
  pass <- achieved <= width + 1e-8
  stopifnot(
    "achieved half-width must not exceed the requested width" = pass
  )

  if (method != "wald" && n > 2) {
    achieved_minus_1 <- independent_max_half_width(n - 1, p, alpha, method)
    stopifnot(
      "n - 1 must fail the precision target (minimality)" = achieved_minus_1 > width - 1e-8
    )
  }

  stopifnot(
    "package achieved_maximum_half_width must match the independent recomputation" =
      isTRUE(all.equal(res$diagnostics$achieved_maximum_half_width, achieved, tolerance = 1e-6))
  )

  results[[i]] <- data.frame(
    p = p, width = width, alpha = alpha, method = method, n = n,
    achieved_max_half_width = achieved, pass = pass,
    expected_events = res$diagnostics$expected_events,
    probability_zero_events = res$diagnostics$probability_zero_events
  )
}

validation_results <- do.call(rbind, results)
message("All invariants passed for ", nrow(validation_results), " (p, width, alpha, method) combinations.")
print(validation_results)

invisible(validation_results)
