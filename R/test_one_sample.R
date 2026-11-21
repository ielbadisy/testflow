#' Test one numeric sample against a reference value
#' @param data A data frame.
#' @param outcome Numeric outcome column.
#' @param mu Reference value.
#' @param alternative Alternative hypothesis.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_one_sample`. The object is
#' a list containing the cleaned data, descriptive statistics, assumption
#' checks, recommended test, primary test result with null hypothesis,
#' alternative test results, effect size, optional `ggplot`, original call, and
#' report text. `$alternative_tests$median_ci` holds the distribution-free
#' confidence interval for the median: with \eqn{k} the \eqn{\alpha/2}
#' quantile of Binomial(\eqn{n}, 0.5), the interval is
#' \eqn{[x_{(k)}, x_{(n-k+1)}]} and its achieved coverage
#' \eqn{1 - 2P(B \le k - 1)} (at least \eqn{1 - \alpha}) is in `conf.level`;
#' the bounds are `NA` when \eqn{n} is too small for that coverage.
#' @references
#' Conover, W. J. (1999). \emph{Practical Nonparametric Statistics} (3rd ed.).
#' Wiley. Section 3.2.
#'

#' Gosset, W. S. (1908). The probable error of a mean. \emph{Biometrika},
#' 6(1), 1-25.
#'
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods.
#' \emph{Biometrics Bulletin}, 1(6), 80-83.
#' @export
testonesample <- function(data, outcome, mu = 0, alternative = c("two.sided", "less", "greater"), alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  alternative <- match.arg(alternative)
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  df <- drop_missing(data, outcome_nm, na.rm = na.rm)
  normality <- check_normality(df, outcome_nm, alpha = alpha)
  recommendation <- recommend_one_sample(normality)
  xval <- df[[outcome_nm]]
  symmetry <- check_symmetry(xval - mu, alpha)
  independence <- check_independence_note()
  ttest <- stats::t.test(xval, mu = mu, alternative = alternative)
  wilcox <- stats::wilcox.test(xval, mu = mu, alternative = alternative, exact = FALSE)
  signs <- stats::binom.test(sum(xval > mu), sum(xval != mu), p = 0.5, alternative = alternative)
  primary <- if (recommendation == "One-sample t-test") ttest else wilcox
  effect <- cohens_d_one_sample(xval, mu)
  plt <- if (plot) make_plot("one_sample", df, outcome_nm, recommended = recommendation, primary = primary, effect = effect, extra = list(mu = mu, median_ci = median_ci_binom(xval, alpha))) else NULL
  h0 <- h0_mean_equal(outcome_nm)
  out <- new_testflow("one_sample", "one numerical sample", outcome_nm, data = df, descriptives = descriptives_numeric(df, outcome_nm), assumptions = assumption_checks(independence, normality, symmetry), recommended = list(test = recommendation, rationale = recommend_rationale(recommendation, normality, alpha = alpha)), primary_test = add_null_hypothesis(safe_tidy_htest(primary, recommendation), h0), alternative_tests = list(t_test = add_null_hypothesis(safe_tidy_htest(ttest, "One-sample t-test"), h0), wilcox = add_null_hypothesis(safe_tidy_htest(wilcox, "Wilcoxon signed-rank test"), h0), sign = add_null_hypothesis(safe_tidy_htest(signs, "Sign test"), paste0("H0: observations are equally likely to be above or below the reference value.")), median_ci = median_ci_binom(xval, alpha)), effect_size = effect, plot = plt, call = match.call(), subclass = "one_sample")
  out$interpretation <- make_report(out, alpha)
  out
}

# Order-statistic (binomial) CI for the median; same as DescTools::MedianCI(method = "exact").
median_ci_binom <- function(x, alpha = 0.05) {
  n <- length(x); s <- sort(x); k <- stats::qbinom(alpha / 2, n, 0.5)
  ok <- k >= 1
  .bt(method = "Order-statistic (binomial) CI for the median", estimate = stats::median(x),
      conf.low = if (ok) s[k] else NA_real_, conf.high = if (ok) s[n - k + 1] else NA_real_,
      conf.level = if (ok) 1 - 2 * stats::pbinom(k - 1, n, 0.5) else NA_real_, rank = k)
}
