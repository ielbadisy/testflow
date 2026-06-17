#' Run a factorial ANOVA workflow
#' @param formula A formula such as `outcome ~ factor1 * factor2`, or a data frame when using data-first style.
#' @param data A data frame, or the outcome column when using data-first style.
#' @param factors Factor columns selected with tidyselect syntax. Optional when using formula style.
#' @param alpha Significance level.
#' @param type ANOVA type placeholder for future car integration.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @references
#' Fisher, R. A. (1925). \emph{Statistical Methods for Research Workers}.
#' Oliver and Boyd.
#'
#' Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral
#' Sciences} (2nd ed.). Lawrence Erlbaum.
#' @export
test_factorial <- function(formula, data, factors = NULL, alpha = 0.05, type = 2, plot = TRUE, na.rm = TRUE) {
  first_expr <- substitute(formula)
  second_expr <- substitute(data)
  data_obj <- resolve_data_first_or_formula(formula, data)
  if (inherits(first_expr, "formula") || (is.call(first_expr) && identical(first_expr[[1]], as.name("~")))) {
    outcome_nm <- all.vars(first_expr[[2]])[1]
    factor_nms <- all.vars(first_expr[[3]])
    formula_obj <- stats::as.formula(deparse(first_expr))
  } else if (inherits(second_expr, "formula") || (is.call(second_expr) && identical(second_expr[[1]], as.name("~")))) {
    outcome_nm <- all.vars(second_expr[[2]])[1]
    factor_nms <- all.vars(second_expr[[3]])
    formula_obj <- stats::as.formula(deparse(second_expr))
  } else {
    outcome_nm <- rlang::as_name(second_expr)
    factor_nms <- tidyselect_names(data_obj, {{ factors }})
    formula_obj <- stats::as.formula(paste(outcome_nm, "~", paste(factor_nms, collapse = " * ")))
  }
  df <- drop_missing(data_obj, c(outcome_nm, factor_nms), na.rm = na.rm)
  df[factor_nms] <- lapply(df[factor_nms], as.factor)
  fit <- stats::aov(formula_obj, data = df)
  tab <- broom::tidy(fit)
  residual_df <- tibble::tibble(.resid = stats::residuals(fit))
  normality <- check_normality(residual_df, ".resid", alpha = alpha)
  levene <- check_variance_homogeneity(df, outcome_nm, factor_nms[1], alpha)
  balanced <- assumption_check("Balanced design", "not required", ifelse(length(unique(table(df[factor_nms]))) > 1, "Cell sizes are unbalanced; the workflow still reports the design.", "Cell sizes are balanced."))
  effect <- eta_squared_aov(fit)
  primary <- tab |> dplyr::filter(.data$term != "Residuals") |> dplyr::slice(1) |> dplyr::transmute(method = "Factorial ANOVA", statistic = .data$statistic, parameter = .data$df, p.value = .data$p.value)
  plt <- if (plot && length(factor_nms) >= 2) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[factor_nms[1]]], y = .data[[outcome_nm]], color = .data[[factor_nms[2]]], group = .data[[factor_nms[2]]])) +
      ggplot2::stat_summary(fun = mean, geom = "line") +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 2) +
      ggplot2::labs(title = "Factorial ANOVA workflow", subtitle = plot_subtitle("Factorial ANOVA", list(p.value = min(primary$p.value, na.rm = TRUE))), x = factor_nms[1], y = outcome_nm, color = factor_nms[2]) +
      ggplot2::theme_minimal()
  } else if (plot) make_plot("groups", df, outcome_nm, factor_nms[1], "Factorial ANOVA", list(p.value = min(primary$p.value, na.rm = TRUE)), effect) else NULL
  h0 <- h0_mean_equal(outcome_nm, paste(factor_nms, collapse = ", "))
  out <- new_testflow("factorial", "factorial design", outcome_nm, paste(factor_nms, collapse = ", "), data = df, descriptives = descriptives_numeric(df, outcome_nm, factor_nms[1]), assumptions = assumption_checks(check_independence_note(), normality, levene, balanced), recommended = list(test = "Factorial ANOVA", rationale = "Primary workflow for factorial numeric outcomes."), primary_test = add_null_hypothesis(primary, h0), alternative_tests = list(anova_table = tab), effect_size = effect, plot = plt, call = match.call(), subclass = "factorial")
  out$interpretation <- make_report(out, alpha)
  out
}
