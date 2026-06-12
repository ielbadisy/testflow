#' Test one numeric sample against a reference value
#' @param data A data frame.
#' @param outcome Numeric outcome column.
#' @param mu Reference value.
#' @param alternative Alternative hypothesis.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_one_sample <- function(data, outcome, mu = 0, alternative = c("two.sided", "less", "greater"), alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  alternative <- match.arg(alternative)
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  df <- drop_missing(data, outcome_nm, na.rm = na.rm)
  normality <- check_normality(df, outcome_nm, alpha = alpha)
  recommendation <- recommend_one_sample(normality)
  xval <- df[[outcome_nm]]
  ttest <- stats::t.test(xval, mu = mu, alternative = alternative)
  wilcox <- stats::wilcox.test(xval, mu = mu, alternative = alternative, exact = FALSE)
  signs <- stats::binom.test(sum(xval > mu), sum(xval != mu), p = 0.5, alternative = alternative)
  primary <- if (recommendation == "One-sample t-test") ttest else wilcox
  effect <- cohens_d_one_sample(xval, mu)
  plt <- if (plot) make_plot("one_sample", df, outcome_nm, recommended = recommendation, primary = primary, effect = effect, extra = list(mu = mu)) else NULL
  h0 <- h0_mean_equal(outcome_nm)
  out <- new_testflow("one_sample", "one numerical sample", outcome_nm, data = df, descriptives = descriptives_numeric(df, outcome_nm), assumptions = list("Normality" = normality), recommended = list(test = recommendation), primary_test = add_null_hypothesis(safe_tidy_htest(primary, recommendation), h0), alternative_tests = list(t_test = add_null_hypothesis(safe_tidy_htest(ttest, "One-sample t-test"), h0), wilcox = add_null_hypothesis(safe_tidy_htest(wilcox, "Wilcoxon signed-rank test"), h0), sign = add_null_hypothesis(safe_tidy_htest(signs, "Sign test"), paste0("H0: observations are equally likely to be above or below the reference value."))), effect_size = effect, plot = plt, call = match.call(), subclass = "one_sample")
  out$interpretation <- make_report(out, alpha)
  out
}
