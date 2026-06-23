#' Compare paired before and after numeric measurements
#' @param formula A formula such as `after ~ before`, or a data frame when using data-first style.
#' @param data A data frame, or the before column when using data-first style.
#' @param after After column. Optional when using formula style.
#' @param alternative Alternative hypothesis.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_paired <- function(formula, data, after = NULL, alternative = c("two.sided", "less", "greater"), alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  alternative <- match.arg(alternative)
  first_expr <- substitute(formula)
  second_expr <- substitute(data)
  is_formula_call <- inherits(first_expr, "formula") ||
    (is.call(first_expr) && identical(first_expr[[1]], as.name("~"))) ||
    inherits(second_expr, "formula") ||
    (is.call(second_expr) && identical(second_expr[[1]], as.name("~")))
  vars <- resolve_formula_pair(first_expr, second_expr, substitute(after), missing(after), rhs_label = "before")
  if (is_formula_call) {
    after_nm <- vars$outcome
    before_nm <- vars$group
  } else {
    before_nm <- vars$outcome
    after_nm <- vars$group
  }
  data_obj <- resolve_data_first_or_formula(formula, data)
  df <- drop_missing(data_obj, c(before_nm, after_nm), na.rm = na.rm)
  df$row_id <- seq_len(nrow(df))
  df$diff <- df[[after_nm]] - df[[before_nm]]
  normality <- check_normality(df, "diff", alpha = alpha)
  symmetry <- assumption_check(
    "Symmetry of paired differences",
    ifelse(normality$status[1] == "acceptable", "not checked", "warning"),
    ifelse(normality$status[1] == "acceptable", "Normality made the symmetry check unnecessary.", "Wilcoxon signed-rank assumes approximate symmetry of paired differences.")
  )
  outliers <- check_outliers(df$diff)
  recommendation <- recommend_paired(normality)
  ttest <- stats::t.test(df[[after_nm]], df[[before_nm]], paired = TRUE, alternative = alternative)
  wilcox <- stats::wilcox.test(df[[after_nm]], df[[before_nm]], paired = TRUE, alternative = alternative, exact = FALSE)
  signs <- stats::binom.test(sum(df$diff > 0), sum(df$diff != 0), p = 0.5, alternative = alternative)
  primary <- if (recommendation == "Paired t-test") ttest else wilcox
  effect <- cohens_d_paired(df[[before_nm]], df[[after_nm]])
  long <- tidyr::pivot_longer(df, dplyr::all_of(c(before_nm, after_nm)), names_to = "time", values_to = "value")
  plt <- if (plot) make_plot("paired", long, paste(after_nm, "-", before_nm), recommended = recommendation, primary = primary, effect = effect) else NULL
  h0 <- paste0("H0: the mean or median paired difference (", after_nm, " - ", before_nm, ") equals 0.")
  out <- new_testflow("paired", "paired measurements", paste(after_nm, "-", before_nm), id = "row_id", data = df, descriptives = descriptives_numeric(df, c(before_nm, after_nm, "diff")), assumptions = assumption_checks(check_independence_note("Paired observations from the same subjects are assumed by design."), normality, symmetry, outliers), recommended = list(test = recommendation), primary_test = add_null_hypothesis(safe_tidy_htest(primary, recommendation), h0), alternative_tests = list(paired_t_test = add_null_hypothesis(safe_tidy_htest(ttest, "Paired t-test"), h0), wilcox = add_null_hypothesis(safe_tidy_htest(wilcox, "Wilcoxon signed-rank test"), h0), sign = add_null_hypothesis(safe_tidy_htest(signs, "Sign test"), "H0: positive and negative paired differences are equally likely.")), effect_size = effect, plot = plt, call = match.call(), subclass = "paired")
  out$interpretation <- make_report(out, alpha)
  out
}
