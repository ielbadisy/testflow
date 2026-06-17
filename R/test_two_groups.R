#' Compare a numeric outcome between two independent groups
#' @param formula A formula such as `outcome ~ group`, or a data frame when using pipe/data-first style.
#' @param data A data frame, or an outcome column when using data-first style.
#' @param group Two-level grouping column. Optional when using formula style.
#' @param alternative Alternative hypothesis.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_two_groups <- function(
  formula,
  data,
  group = NULL,
  alternative = c("two.sided", "less", "greater"),
  alpha = 0.05,
  plot = TRUE,
  na.rm = TRUE
) {
  alternative <- match.arg(alternative)
  vars <- resolve_formula_pair(substitute(formula), substitute(data), substitute(group), missing(group))
  outcome_nm <- vars$outcome
  group_nm <- vars$group
  data_obj <- resolve_data_first_or_formula(formula, data)
  df <- drop_missing(data_obj, c(outcome_nm, group_nm), na.rm = na.rm)
  df[[group_nm]] <- as.factor(df[[group_nm]])
  assert_two_groups(df, group_nm)

  normality <- check_normality(df, outcome_nm, group_nm, alpha)
  levene <- check_variance_homogeneity(df, outcome_nm, group_nm, alpha)
  f_test <- check_variance_two_groups(df, outcome_nm, group_nm, alpha)
  outliers <- check_outliers(df[[outcome_nm]])
  independence <- check_independence_note()
  recommendation <- recommend_two_groups(normality, levene)
  formula <- stats::as.formula(paste(outcome_nm, "~", group_nm))

  student <- stats::t.test(formula, data = df, var.equal = TRUE, alternative = alternative)
  welch <- stats::t.test(formula, data = df, var.equal = FALSE, alternative = alternative)
  wilcox <- stats::wilcox.test(formula, data = df, alternative = alternative, exact = FALSE)
  primary <- switch(
    recommendation,
    "Student independent t-test" = student,
    "Welch t-test" = welch,
    "Wilcoxon rank-sum test" = wilcox
  )
  effect <- if (recommendation == "Wilcoxon rank-sum test") rank_biserial_two_groups(df, outcome_nm, group_nm) else cohens_d_independent(df, outcome_nm, group_nm)
  primary_tidy <- add_null_hypothesis(safe_tidy_htest(primary, recommendation), h0_mean_equal(outcome_nm, group_nm))
  plt <- if (plot) make_plot("two_groups", df, outcome_nm, group_nm, recommendation, primary, effect) else NULL

  x <- new_testflow(
    workflow = "two_groups",
    design = "two independent groups",
    outcome = outcome_nm,
    group = group_nm,
    data = df,
    descriptives = descriptives_numeric(df, outcome_nm, group_nm),
    assumptions = assumption_checks(independence, normality, levene, outliers, tibble::add_column(f_test, name = "Variance ratio check", message = "Variance comparison via F test.", .before = 1)),
    recommended = list(test = recommendation, rationale = "Selected from normality and variance assumptions."),
    primary_test = primary_tidy,
    alternative_tests = list(
      student_t_test = add_null_hypothesis(safe_tidy_htest(student, "Student independent t-test"), h0_mean_equal(outcome_nm, group_nm)),
      welch_t_test = add_null_hypothesis(safe_tidy_htest(welch, "Welch t-test"), h0_mean_equal(outcome_nm, group_nm)),
      wilcox = add_null_hypothesis(safe_tidy_htest(wilcox, "Wilcoxon rank-sum test"), h0_mean_equal(outcome_nm, group_nm))
    ),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "two_groups"
  )
  x$interpretation <- make_report(x, alpha)
  x
}
