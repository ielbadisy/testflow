#' Bland-Altman method-comparison workflow
#'
#' Quantifies agreement between two measurement methods applied to the same
#' subjects: the mean difference (bias), the 95% limits of agreement, and
#' confidence intervals for all three.
#'
#' @param data A data frame.
#' @param method1 Column with the first method's measurements.
#' @param method2 Column with the second method's measurements. Differences
#'   are `method1 - method2`.
#' @param alpha Significance level. The limits of agreement are
#'   `bias +/- z(1 - alpha/2) * SD` (1.96 SD when `alpha = 0.05`) and every
#'   CI has level `1 - alpha`.
#' @param plot Logical; include a ggplot object (difference against mean).
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_blandaltman`. The
#' primary result is the one-sample t-test of zero bias (with the bias CI).
#' `$limits` (also `$alternative_tests$limits_table`) holds bias, lower and
#' upper limits of agreement with their CIs. `$alternative_tests` also holds
#' a proportional-bias regression of the difference on the mean.
#' @details
#' With differences \eqn{d_i}, mean \eqn{\bar d}, SD \eqn{s} and
#' \eqn{z = z_{1-\alpha/2}}:
#' \deqn{LoA = \bar d \pm z s}
#' The bias CI is \eqn{\bar d \pm t_{n-1} s/\sqrt n}. The CI of each limit
#' uses the exact-variance form of Bland and Altman (1999),
#' \eqn{SE = s\sqrt{1/n + z^2/(2(n-1))}}, with \eqn{t_{n-1}} quantiles. The
#' 1986 paper (and the BlandAltmanLeh package) uses the approximation
#' \eqn{SE \approx s\sqrt{3/n}}, which is slightly larger for small n.
#'
#' The limits are only meaningful if the differences are roughly normal and
#' do not depend on the magnitude of the measurement. Both are checked
#' (Shapiro-Wilk; slope of difference on mean), and nothing is chosen
#' automatically.
#' @references
#' Bland JM, Altman DG. Statistical methods for assessing agreement between
#' two methods of clinical measurement. *Lancet*. 1986;327(8476):307-310.
#'
#' Bland JM, Altman DG. Measuring agreement in method comparison studies.
#' *Statistical Methods in Medical Research*. 1999;8(2):135-160.
#' @examples
#' set.seed(1)
#' x <- rnorm(40, 120, 15)
#' d <- data.frame(device = x + rnorm(40, 2, 5), reference = x + rnorm(40, 0, 5))
#' testblandaltman(d, device, reference, plot = FALSE)
#' @export
testblandaltman <- function(data, method1, method2, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  m1 <- rlang::as_name(rlang::ensym(method1))
  m2 <- rlang::as_name(rlang::ensym(method2))
  df <- drop_missing(data, c(m1, m2), na.rm = na.rm)
  if (!is.numeric(df[[m1]]) || !is.numeric(df[[m2]])) {
    stop("Both methods must be numeric measurements.", call. = FALSE)
  }
  n <- nrow(df)
  if (n < 3) stop("At least 3 complete pairs are needed.", call. = FALSE)

  df$.diff <- df[[m1]] - df[[m2]]
  df$.mean <- (df[[m1]] + df[[m2]]) / 2
  bias <- mean(df$.diff)
  s <- stats::sd(df$.diff)
  z <- stats::qnorm(1 - alpha / 2)
  tq <- stats::qt(1 - alpha / 2, n - 1)
  se_bias <- s / sqrt(n)
  se_loa <- s * sqrt(1 / n + z^2 / (2 * (n - 1)))
  lower <- bias - z * s
  upper <- bias + z * s

  limits <- .bt(
    quantity = c("Bias", "Lower limit of agreement", "Upper limit of agreement"),
    estimate = c(bias, lower, upper),
    se = c(se_bias, se_loa, se_loa),
    conf.low = c(bias, lower, upper) - tq * c(se_bias, se_loa, se_loa),
    conf.high = c(bias, lower, upper) + tq * c(se_bias, se_loa, se_loa)
  )

  ttest <- stats::t.test(df[[m1]], df[[m2]], paired = TRUE, conf.level = 1 - alpha)
  h0 <- paste0("H0: the mean difference (", m1, " - ", m2, ") equals 0 (no bias).")
  primary <- add_null_hypothesis(safe_tidy_htest(ttest, "Bland-Altman bias (paired t-test)"), h0)

  fit <- stats::lm(.diff ~ .mean, data = df)
  slope <- summary(fit)$coefficients[2, , drop = FALSE]
  prop <- .bt(
    term = "slope of difference on mean",
    estimate = slope[1, 1], std.error = slope[1, 2], statistic = slope[1, 3], p.value = slope[1, 4],
    null_hypothesis = "H0: the difference does not change with the magnitude of measurement."
  )
  normality <- check_normality(df, ".diff", alpha = alpha)
  normality$name <- "Normality of differences"
  prop_check <- assumption_check(
    "No proportional bias",
    if (slope[1, 4] < alpha) "not acceptable" else "acceptable",
    if (slope[1, 4] < alpha) "The difference changes with the mean; constant limits may mislead (consider log scale or regression-based limits)." else "No evidence that the difference depends on the magnitude.",
    method = "Regression of difference on mean", statistic = slope[1, 3], p_value = slope[1, 4]
  )

  effect <- .bt(name = "Bias (mean difference)", estimate = bias, magnitude = NA_character_)

  plt <- if (plot) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data$.mean, y = .data$.diff)) +
      ggplot2::geom_hline(yintercept = c(lower, upper), linetype = "dashed", colour = "#B23A48") +
      ggplot2::geom_hline(yintercept = bias, colour = "#1F2933") +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
      ggplot2::geom_point(alpha = 0.8, colour = "#3E6A9E") +
      ggplot2::labs(
        title = "Bland-Altman plot",
        subtitle = paste0("Bias ", format_stat(bias), "; limits of agreement [", format_stat(lower), ", ", format_stat(upper), "]"),
        x = paste0("Mean of ", m1, " and ", m2), y = paste0(m1, " - ", m2)
      ) +
      ggplot2::theme_minimal()
  } else NULL

  out <- new_testflow(
    "blandaltman",
    "method comparison (Bland-Altman)",
    m1,
    m2,
    data = df,
    descriptives = descriptives_numeric(df, c(m1, m2, ".diff")),
    assumptions = assumption_checks(check_independence_note("Each subject is measured once by each method."), normality, prop_check),
    recommended = list(test = "Bland-Altman limits of agreement", rationale = "Two methods measuring the same quantity on the same subjects; agreement is judged against a clinically acceptable difference, not by a test."),
    primary_test = primary,
    alternative_tests = list(limits_table = limits, proportional_bias = prop),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "blandaltman"
  )
  out$limits <- limits
  out$interpretation <- paste0(
    make_report(out, alpha),
    " The ", format(100 * (1 - alpha)), "% limits of agreement were [", format_stat(lower), ", ", format_stat(upper), "]",
    " (lower CI ", format_stat(limits$conf.low[2]), " to ", format_stat(limits$conf.high[2]),
    "; upper CI ", format_stat(limits$conf.low[3]), " to ", format_stat(limits$conf.high[3]), ")."
  )
  out
}
