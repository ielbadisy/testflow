#' Compare a numeric outcome across more than two groups
#' @param formula A formula such as `outcome ~ group`, or a data frame when using pipe/data-first style.
#' @param data A data frame, or an outcome column when using data-first style.
#' @param group Grouping column. Optional when using formula style.
#' @param alpha Significance level.
#' @param posthoc Logical; compute post-hoc comparisons.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_groups <- function(formula, data, group = NULL, alpha = 0.05, posthoc = TRUE, plot = TRUE, na.rm = TRUE) {
  vars <- resolve_formula_pair(substitute(formula), substitute(data), substitute(group), missing(group))
  outcome_nm <- vars$outcome
  group_nm <- vars$group
  data_obj <- resolve_data_first_or_formula(formula, data)
  df <- drop_missing(data_obj, c(outcome_nm, group_nm), na.rm = na.rm)
  df[[group_nm]] <- as.factor(df[[group_nm]])
  if (dplyr::n_distinct(df[[group_nm]]) < 3) stop("`group` must contain at least three groups.", call. = FALSE)
  normality <- check_normality(df, outcome_nm, group_nm, alpha)
  levene <- check_variance_homogeneity(df, outcome_nm, group_nm, alpha)
  bartlett <- check_bartlett(df, outcome_nm, group_nm, alpha)
  independence <- check_independence_note()
  outliers <- check_outliers(df[[outcome_nm]])
  recommendation <- recommend_groups(normality, levene)
  formula <- stats::as.formula(paste(outcome_nm, "~", group_nm))
  aov_fit <- stats::aov(formula, data = df)
  aov_test <- broom::tidy(aov_fit) |> dplyr::slice(1) |> dplyr::rename(p.value = "p.value")
  welch <- stats::oneway.test(formula, data = df, var.equal = FALSE)
  kruskal <- stats::kruskal.test(formula, data = df)
  primary <- switch(recommendation, "One-way ANOVA" = list_obj_from_tidy(aov_test, "One-way ANOVA"), "Welch ANOVA" = welch, "Kruskal-Wallis test" = kruskal)
  h0 <- h0_mean_equal(outcome_nm, group_nm)
  primary_tidy <- if (inherits(primary, "htest")) safe_tidy_htest(primary, recommendation) else primary
  primary_tidy <- add_null_hypothesis(primary_tidy, h0)
  ph <- if (posthoc) posthoc_groups(df, outcome_nm, group_nm, recommendation, alpha = alpha) else NULL
  effect <- if (recommendation == "Kruskal-Wallis test") {
    h <- unname(kruskal$statistic); n <- nrow(df); k <- dplyr::n_distinct(df[[group_nm]])
    tibble::tibble(name = "Kruskal epsilon squared", estimate = (h - k + 1) / (n - k), magnitude = magnitude_eta2((h - k + 1) / (n - k)))
  } else eta_squared_aov(aov_fit) |> dplyr::slice(1)
  plt <- if (plot) make_plot("groups", df, outcome_nm, group_nm, recommendation, if (inherits(primary, "htest")) primary else list(p.value = primary_tidy$p.value[1]), effect) else NULL
  out <- new_testflow("groups", "more than two independent groups", outcome_nm, group_nm, data = df, descriptives = descriptives_numeric(df, outcome_nm, group_nm), assumptions = assumption_checks(independence, normality, levene, bartlett, outliers), recommended = list(test = recommendation), primary_test = primary_tidy, alternative_tests = list(anova = add_null_hypothesis(aov_test, h0), welch = add_null_hypothesis(safe_tidy_htest(welch, "Welch ANOVA"), h0), kruskal = add_null_hypothesis(safe_tidy_htest(kruskal, "Kruskal-Wallis test"), h0)), posthoc = ph, effect_size = effect, plot = plt, call = match.call(), subclass = "groups")
  out$interpretation <- make_report(out, alpha)
  out
}

list_obj_from_tidy <- function(x, method) {
  tibble::tibble(method = method, statistic = x$statistic[1], parameter = x$df[1], p.value = x$p.value[1])
}
