#' Run a factorial ANOVA workflow
#' @param formula A formula such as `outcome ~ factor1 * factor2`, or a data frame when using data-first style.
#' @param data A data frame, or the outcome column when using data-first style.
#' @param factors Factor columns selected with tidyselect syntax. Optional when using formula style.
#' @param alpha Significance level.
#' @param type Sums-of-squares type: `1` (sequential, base `aov()`), `2`, or
#'   `3` (via `car::Anova()`, required for `type = 2`/`3`). For unbalanced
#'   designs these can give materially different p-values; `type = 2` is a
#'   reasonable default when there is no strong prior reason to test factors
#'   in a particular order. `type = 3` requires sum-to-zero contrasts, which
#'   this function sets automatically.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_factorial`. The object is a
#' list containing the cleaned data, descriptive statistics, residual and
#' variance assumption checks, recommended factorial ANOVA, primary ANOVA term
#' result with null hypothesis, ANOVA table, effect size, optional `ggplot`,
#' original call, and report text.
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
  warn_if_two_groups_for_factorial(factor_nms)
  stopifnot(type %in% c(1, 2, 3))
  df <- drop_missing(data_obj, c(outcome_nm, factor_nms), na.rm = na.rm)
  df[factor_nms] <- lapply(df[factor_nms], as.factor)
  fit <- stats::aov(formula_obj, data = df)
  ss_tab <- if (identical(type, 1)) {
    stats::anova(fit)
  } else {
    if (!requireNamespace("car", quietly = TRUE)) {
      stop("`type = 2` and `type = 3` sums of squares require the 'car' package. Install it, or use `type = 1`.", call. = FALSE)
    }
    if (identical(type, 3)) {
      # Type III SS are only interpretable with sum-to-zero contrasts;
      # under the default treatment contrasts they test a different (usually
      # meaningless) hypothesis about a reference-level intercept.
      old_contrasts <- options(contrasts = c("contr.sum", "contr.poly"))
      on.exit(options(old_contrasts), add = TRUE)
      lm_fit <- stats::lm(formula_obj, data = df)
    } else {
      lm_fit <- stats::lm(formula_obj, data = df)
    }
    car::Anova(lm_fit, type = type)
  }
  tab <- broom::tidy(ss_tab)
  tab <- tab[tab$term != "(Intercept)", , drop = FALSE]
  residual_df <- tibble::tibble(.resid = stats::residuals(fit))
  residual_normality_test <- stats::shapiro.test(residual_df$.resid)
  normality <- assumption_check(
    "Normality of residuals",
    ifelse(residual_normality_test$p.value >= alpha, "acceptable", "not acceptable"),
    ifelse(residual_normality_test$p.value >= alpha, "Residuals appear approximately normal.", "Residuals deviate from normality."),
    method = "Shapiro-Wilk",
    statistic = unname(residual_normality_test$statistic),
    p_value = residual_normality_test$p.value
  )
  levene_test <- check_variance_homogeneity(df, outcome_nm, factor_nms[1], alpha)
  levene <- assumption_check(
    "Variance homogeneity",
    levene_test$status[1],
    ifelse(levene_test$status[1] == "acceptable", "Variance homogeneity looks reasonable.", "Variance homogeneity may be violated."),
    method = levene_test$method[1],
    statistic = levene_test$statistic[1],
    p_value = levene_test$p[1],
    details = paste0("Df1=", levene_test$df1[1], "; Df2=", levene_test$df2[1])
  )
  type_label <- c(`1` = "Type I (sequential)", `2` = "Type II", `3` = "Type III")[[as.character(type)]]
  balanced <- assumption_check("Balanced design", "not required", ifelse(length(unique(table(df[factor_nms]))) > 1, paste0("Cell sizes are unbalanced; ", type_label, " sums of squares are used so the terms are not order-dependent.", if (identical(type, 1)) " Type I sums of squares depend on the order factors are listed in the formula for unbalanced designs; consider `type = 2` or `type = 3`." else ""), "Cell sizes are balanced; Type I, II, and III sums of squares agree."))
  effect <- eta_squared_aov(fit, tab = ss_tab)
  primary <- tab |> dplyr::filter(.data$term != "Residuals") |> dplyr::slice(1) |> dplyr::transmute(method = paste0("Factorial ANOVA (", type_label, ")"), statistic = .data$statistic, parameter = .data$df, p.value = .data$p.value)
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
