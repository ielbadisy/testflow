#' Run a factorial ANOVA workflow
#' @param data A data frame.
#' @param outcome Numeric outcome column.
#' @param factors Factor columns selected with tidyselect syntax.
#' @param alpha Significance level.
#' @param type ANOVA type placeholder for future car integration.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_factorial <- function(data, outcome, factors, alpha = 0.05, type = 2, plot = TRUE, na.rm = TRUE) {
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  factor_nms <- tidyselect_names(data, {{ factors }})
  df <- drop_missing(data, c(outcome_nm, factor_nms), na.rm = na.rm)
  df[factor_nms] <- lapply(df[factor_nms], as.factor)
  formula <- stats::as.formula(paste(outcome_nm, "~", paste(factor_nms, collapse = " * ")))
  fit <- stats::aov(formula, data = df)
  tab <- broom::tidy(fit)
  residual_df <- tibble::tibble(.resid = stats::residuals(fit))
  normality <- check_normality(residual_df, ".resid", alpha = alpha)
  levene <- check_variance_homogeneity(df, outcome_nm, factor_nms[1], alpha)
  effect <- eta_squared_aov(fit)
  primary <- tab |> dplyr::filter(.data$term != "Residuals") |> dplyr::slice(1) |> dplyr::transmute(method = "Factorial ANOVA", statistic = .data$statistic, parameter = .data$df, p.value = .data$p.value)
  plt <- if (plot && length(factor_nms) >= 2) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[factor_nms[1]]], y = .data[[outcome_nm]], color = .data[[factor_nms[2]]], group = .data[[factor_nms[2]]])) +
      ggplot2::stat_summary(fun = mean, geom = "line") +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 2) +
      ggplot2::labs(title = "Factorial ANOVA workflow", subtitle = plot_subtitle("Factorial ANOVA", list(p.value = min(primary$p.value, na.rm = TRUE))), x = factor_nms[1], y = outcome_nm, color = factor_nms[2]) +
      ggplot2::theme_minimal()
  } else if (plot) make_plot("groups", df, outcome_nm, factor_nms[1], "Factorial ANOVA", list(p.value = min(primary$p.value, na.rm = TRUE)), effect) else NULL
  out <- new_testflow("factorial", "factorial design", outcome_nm, paste(factor_nms, collapse = ", "), data = df, descriptives = descriptives_numeric(df, outcome_nm, factor_nms[1]), assumptions = list("Residual normality" = normality, "Homogeneity of variance" = levene), recommended = list(test = "Factorial ANOVA", rationale = "Primary workflow for factorial numeric outcomes."), primary_test = primary, alternative_tests = list(anova_table = tab), effect_size = effect, plot = plt, call = match.call(), subclass = "factorial")
  out$interpretation <- make_report(out, alpha)
  out
}
