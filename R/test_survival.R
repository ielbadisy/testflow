#' Kaplan-Meier and log-rank survival workflow
#'
#' @param formula A formula such as `Surv(time, status) ~ group`, with a
#'   two-level grouping factor.
#' @param data A data frame.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_survival`. The object is
#' a list containing the cleaned data, per-group descriptive survival
#' statistics, the log-rank test as the primary result with null hypothesis,
#' the Kaplan-Meier curve data, a companion univariate hazard ratio as the
#' effect size, optional `ggplot`, original call, and report text.
#' @details
#' For groups A and B, the log-rank test compares observed to expected
#' events under the null of equal survival distributions:
#' \deqn{\chi^2 = \sum_{g \in \{A,B\}} \frac{(O_g - E_g)^2}{V_g} \sim \chi^2_1}
#'
#' The companion effect size is the hazard ratio from a univariate Cox model
#' on the same grouping factor, \eqn{HR = e^{\hat\beta}}, reported with its
#' Wald confidence interval; the hazard ratio is not itself part of the
#' log-rank test and can disagree with it under strongly non-proportional
#' hazards.
#' @references
#' Kaplan EL, Meier P. Nonparametric estimation from incomplete
#' observations. *Journal of the American Statistical Association*.
#' 1958;53(282):457-481.
#'
#' Mantel N. Evaluation of survival data and two new rank order statistics
#' arising in its consideration. *Cancer Chemotherapy Reports*.
#' 1966;50(3):163-170.
#' @export
test_survival <- function(formula, data, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  formula_obj <- stats::as.formula(deparse(formula))
  lhs <- formula_obj[[2]]
  time_nm <- deparse(lhs[[2]])
  status_nm <- deparse(lhs[[3]])
  group_nm <- all.vars(formula_obj[[3]])[1]

  df <- drop_missing(data, c(time_nm, status_nm, group_nm), na.rm = na.rm)
  assert_two_groups(df, group_nm)
  df[[group_nm]] <- droplevels(as.factor(df[[group_nm]]))

  fit_km <- survival::survfit(formula_obj, data = df)
  fit_logrank <- survival::survdiff(formula_obj, data = df)
  fit_cox <- survival::coxph(formula_obj, data = df)

  km_table <- broom::tidy(fit_km)
  km_table$strata <- sub(paste0("^", group_nm, "="), "", km_table$strata)
  summary_table <- summary(fit_km)$table
  descriptives <- tibble::tibble(
    group = sub(paste0("^", group_nm, "="), "", rownames(summary_table)),
    n = unname(summary_table[, "records"]),
    events = unname(summary_table[, "events"]),
    median_survival = unname(summary_table[, "median"])
  )

  chisq <- fit_logrank$chisq
  lr_df <- length(fit_logrank$n) - 1
  lr_p <- stats::pchisq(chisq, lr_df, lower.tail = FALSE)

  hr <- unname(exp(stats::coef(fit_cox))[1])
  hr_ci <- exp(stats::confint(fit_cox))[1, ]
  effect <- tibble::tibble(name = "Hazard ratio", estimate = hr, magnitude = NA_character_)

  primary <- tibble::tibble(
    method = "Log-rank test",
    statistic = chisq,
    parameter = lr_df,
    p.value = lr_p,
    conf.low = unname(hr_ci[1]),
    conf.high = unname(hr_ci[2])
  )
  h0 <- paste0("H0: the survival distributions are equal across levels of ", group_nm, ".")

  plt <- if (plot) {
    ggplot2::ggplot(km_table, ggplot2::aes(x = .data$time, y = .data$estimate, color = .data$strata)) +
      ggplot2::geom_step(linewidth = 0.8) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high, fill = .data$strata), alpha = 0.12, color = NA) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(
        title = "Survival workflow",
        subtitle = paste0("Log-rank test, p = ", format_p(lr_p), "; hazard ratio = ", format_stat(hr)),
        x = time_nm, y = "Survival probability", color = group_nm, fill = group_nm
      ) +
      ggplot2::theme_minimal()
  } else NULL

  out <- new_testflow(
    "survival",
    "Kaplan-Meier and log-rank",
    time_nm,
    group_nm,
    data = df,
    descriptives = descriptives,
    assumptions = assumption_checks(check_independence_note("Independent censoring is assumed.")),
    recommended = list(test = "Log-rank test", rationale = "Two-group comparison of survival distributions."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = list(kaplan_meier = km_table, cox = broom::tidy(fit_cox, conf.int = TRUE, exponentiate = TRUE)),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "survival"
  )
  out$interpretation <- make_report(out, alpha)
  out
}

#' Cox proportional hazards regression workflow
#'
#' @param formula A formula such as `Surv(time, status) ~ x1 + x2`.
#' @param data A data frame.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_cox`. The object is a
#' list containing the cleaned data, descriptive statistics, a proportional-
#' hazards assumption check, the overall likelihood-ratio test as the
#' primary result with null hypothesis, the per-term hazard-ratio
#' coefficient table, the concordance index as the effect size, optional
#' `ggplot`, original call, and report text.
#' @details
#' With predictors \eqn{x_1, \ldots, x_p}, the Cox model leaves the baseline
#' hazard \eqn{h_0(t)} unspecified:
#' \deqn{h(t \mid x) = h_0(t)\exp(\beta_1 x_1 + \cdots + \beta_p x_p)}
#'
#' Coefficients are reported as hazard ratios \eqn{HR_j = e^{\beta_j}}. The
#' overall model test is the likelihood-ratio test against the model with no
#' predictors, on \eqn{p} degrees of freedom (or more, for multi-level
#' factor terms).
#'
#' The proportional-hazards assumption is checked via the correlation
#' between scaled Schoenfeld residuals and time (Grambsch & Therneau, 1994);
#' a significant test suggests a predictor's effect on the hazard changes
#' over time.
#' @references
#' Cox DR. Regression models and life-tables. *Journal of the Royal
#' Statistical Society: Series B*. 1972;34(2):187-202.
#'
#' Grambsch PM, Therneau TM. Proportional hazards tests and diagnostics
#' based on weighted residuals. *Biometrika*. 1994;81(3):515-526.
#' @export
test_cox <- function(formula, data, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  formula_obj <- stats::as.formula(deparse(formula))
  lhs <- formula_obj[[2]]
  time_nm <- deparse(lhs[[2]])
  status_nm <- deparse(lhs[[3]])
  predictor_nms <- all.vars(formula_obj[[3]])

  df <- drop_missing(data, c(time_nm, status_nm, predictor_nms), na.rm = na.rm)
  fit <- survival::coxph(formula_obj, data = df)
  g <- broom::glance(fit)
  logtest <- summary(fit)$logtest

  coefs_hr <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  names(coefs_hr)[names(coefs_hr) == "estimate"] <- "hazard_ratio"
  coefs_log <- broom::tidy(fit, conf.int = TRUE)

  zph <- survival::cox.zph(fit)
  zph_global <- zph$table["GLOBAL", ]
  ph_check <- assumption_check(
    "Proportional hazards",
    ifelse(zph_global["p"] >= alpha, "acceptable", "not acceptable"),
    ifelse(zph_global["p"] >= alpha, "No evidence against proportional hazards.", "Scaled Schoenfeld residuals are associated with time; proportional hazards may be violated."),
    method = "Schoenfeld residual test (Grambsch-Therneau)",
    statistic = unname(zph_global["chisq"]),
    p_value = unname(zph_global["p"])
  )

  effect <- tibble::tibble(name = "Concordance index", estimate = unname(g$concordance), magnitude = NA_character_)
  primary <- tibble::tibble(
    method = "Cox proportional hazards regression (likelihood-ratio test)",
    statistic = unname(logtest["test"]),
    parameter = unname(logtest["df"]),
    p.value = unname(logtest["pvalue"])
  )
  h0 <- paste0("H0: none of ", paste(predictor_nms, collapse = ", "), " are associated with the hazard (all log hazard ratios equal 0).")

  plt <- if (plot) {
    hr_df <- coefs_hr
    ggplot2::ggplot(hr_df, ggplot2::aes(x = .data$hazard_ratio, y = stats::reorder(.data$term, .data$hazard_ratio))) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#8A8F98") +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = .data$conf.low, xmax = .data$conf.high), color = "#4C78A8") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Cox regression workflow", subtitle = "Hazard ratios with 95% confidence intervals", x = "Hazard ratio (log scale)", y = NULL) +
      ggplot2::theme_minimal()
  } else NULL

  out <- new_testflow(
    "cox",
    "Cox proportional hazards regression",
    time_nm,
    paste(predictor_nms, collapse = ", "),
    data = df,
    descriptives = descriptives_numeric(df, intersect(predictor_nms, names(df)[vapply(df, is.numeric, logical(1))])),
    assumptions = assumption_checks(check_independence_note("Independent censoring is assumed."), ph_check),
    recommended = list(test = "Cox proportional hazards regression", rationale = "Time-to-event outcome modeled as a function of the predictors."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = list(coefficients = coefs_log, hazard_ratios = coefs_hr),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "cox"
  )
  out$interpretation <- make_report(out, alpha)
  out
}
