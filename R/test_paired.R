#' Compare paired before and after numeric measurements
#' @param data A data frame.
#' @param before Before column.
#' @param after After column.
#' @param alternative Alternative hypothesis.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_paired <- function(data, before, after, alternative = c("two.sided", "less", "greater"), alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  alternative <- match.arg(alternative)
  before_nm <- rlang::as_name(rlang::ensym(before))
  after_nm <- rlang::as_name(rlang::ensym(after))
  df <- drop_missing(data, c(before_nm, after_nm), na.rm = na.rm)
  df$row_id <- seq_len(nrow(df))
  df$diff <- df[[after_nm]] - df[[before_nm]]
  normality <- check_normality(df, "diff", alpha = alpha)
  recommendation <- recommend_paired(normality)
  ttest <- stats::t.test(df[[after_nm]], df[[before_nm]], paired = TRUE, alternative = alternative)
  wilcox <- stats::wilcox.test(df[[after_nm]], df[[before_nm]], paired = TRUE, alternative = alternative, exact = FALSE)
  signs <- stats::binom.test(sum(df$diff > 0), sum(df$diff != 0), p = 0.5, alternative = alternative)
  primary <- if (recommendation == "Paired t-test") ttest else wilcox
  effect <- cohens_d_paired(df[[before_nm]], df[[after_nm]])
  long <- tidyr::pivot_longer(df, dplyr::all_of(c(before_nm, after_nm)), names_to = "time", values_to = "value")
  plt <- if (plot) make_plot("paired", long, paste(after_nm, "-", before_nm), recommended = recommendation, primary = primary, effect = effect) else NULL
  out <- new_testflow("paired", "paired measurements", paste(after_nm, "-", before_nm), id = "row_id", data = df, descriptives = descriptives_numeric(df, c(before_nm, after_nm, "diff")), assumptions = list("Normality of paired differences" = normality), recommended = list(test = recommendation), primary_test = safe_tidy_htest(primary, recommendation), alternative_tests = list(paired_t_test = safe_tidy_htest(ttest, "Paired t-test"), wilcox = safe_tidy_htest(wilcox, "Wilcoxon signed-rank test"), sign = safe_tidy_htest(signs, "Sign test")), effect_size = effect, plot = plt, call = match.call(), subclass = "paired")
  out$interpretation <- make_report(out, alpha)
  out
}
