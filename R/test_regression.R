#' Multiple linear regression workflow
#'
#' @param formula A formula such as `outcome ~ x1 + x2 + x3`.
#' @param data A data frame.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_linear_regression`. The
#' object is a list containing the cleaned data, descriptive statistics,
#' residual/homoscedasticity/multicollinearity assumption checks, the overall
#' model F-test as the primary result with null hypothesis, the per-term
#' coefficient table, R squared and adjusted R squared as the effect size,
#' optional `ggplot`, original call, and report text.
#' @details
#' For predictors \eqn{x_1, \ldots, x_p} and \eqn{n} observations, ordinary
#' least squares fits:
#' \deqn{y_i = \beta_0 + \beta_1 x_{1i} + \cdots + \beta_p x_{pi} + \varepsilon_i}
#'
#' The overall model test compares the fitted model to an intercept-only
#' model:
#' \deqn{F = \frac{(SS_{total}-SS_{res})/p}{SS_{res}/(n-p-1)} \sim F_{p,\,n-p-1}}
#'
#' \deqn{R^2 = 1 - \frac{SS_{res}}{SS_{total}}, \qquad
#' R^2_{adj} = 1 - (1-R^2)\frac{n-1}{n-p-1}}
#' @references
#' Draper NR, Smith H. *Applied Regression Analysis* (3rd ed.). Wiley; 1998.
#' @export
test_linear_regression <- function(formula, data, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  formula_obj <- stats::as.formula(deparse(formula))
  outcome_nm <- all.vars(formula_obj[[2]])[1]
  predictor_nms <- all.vars(formula_obj[[3]])
  df <- drop_missing(data, c(outcome_nm, predictor_nms), na.rm = na.rm)
  fit <- stats::lm(formula_obj, data = df)
  g <- broom::glance(fit)
  coefs <- broom::tidy(fit, conf.int = TRUE)

  resid_normality_test <- stats::shapiro.test(stats::residuals(fit))
  normality <- assumption_check(
    "Normality of residuals",
    ifelse(resid_normality_test$p.value >= alpha, "acceptable", "not acceptable"),
    ifelse(resid_normality_test$p.value >= alpha, "Residuals appear approximately normal.", "Residuals deviate from normality."),
    method = "Shapiro-Wilk",
    statistic = unname(resid_normality_test$statistic),
    p_value = resid_normality_test$p.value
  )

  homoscedasticity <- if (requireNamespace("car", quietly = TRUE)) {
    nv <- car::ncvTest(fit)
    assumption_check(
      "Homoscedasticity",
      ifelse(nv$p >= alpha, "acceptable", "not acceptable"),
      ifelse(nv$p >= alpha, "Residual variance looks constant across fitted values.", "Residual variance may depend on the fitted values (heteroscedasticity)."),
      method = "Breusch-Pagan (score test)",
      statistic = unname(nv$ChiSquare),
      p_value = nv$p
    )
  } else {
    assumption_check("Homoscedasticity", "not checked", "Install the 'car' package to run the Breusch-Pagan test.")
  }

  multicollinearity <- if (length(predictor_nms) < 2) {
    assumption_check("Multicollinearity", "not applicable", "Only one predictor; multicollinearity does not apply.")
  } else if (requireNamespace("car", quietly = TRUE)) {
    vifs <- car::vif(fit)
    max_vif <- max(vifs)
    assumption_check(
      "Multicollinearity",
      ifelse(max_vif < 5, "acceptable", "not acceptable"),
      ifelse(max_vif < 5, "Variance inflation factors are all below 5.", "At least one predictor has a variance inflation factor of 5 or higher, indicating multicollinearity."),
      method = "Variance inflation factor",
      statistic = max_vif,
      details = paste0("Max VIF (", names(vifs)[which.max(vifs)], ") = ", format_stat(max_vif))
    )
  } else {
    assumption_check("Multicollinearity", "not checked", "Install the 'car' package to compute variance inflation factors.")
  }

  effect <- dplyr::bind_rows(r_squared_lm(fit), adjusted_r_squared_lm(fit))
  primary <- tibble::tibble(
    method = "Linear regression (overall F-test)",
    statistic = unname(g$statistic),
    parameter = unname(g$df),
    df.error = unname(g$df.residual),
    p.value = unname(g$p.value)
  )
  h0 <- paste0("H0: none of ", paste(predictor_nms, collapse = ", "), " explain variation in ", outcome_nm, " (all slopes equal 0).")

  plt <- if (plot) {
    plot_df <- tibble::tibble(.fitted = stats::fitted(fit), .resid = stats::residuals(fit))
    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$.fitted, y = .data$.resid)) +
      ggplot2::geom_point(alpha = 0.6, color = "#4C78A8") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#F58518") +
      ggplot2::geom_smooth(se = FALSE, method = "loess", formula = y ~ x, color = "#1F2933", linewidth = 0.7) +
      ggplot2::labs(
        title = "Linear regression workflow",
        subtitle = paste0("Residuals vs fitted, R-squared = ", format_stat(g$r.squared)),
        x = "Fitted values", y = "Residuals"
      ) +
      ggplot2::theme_minimal()
  } else NULL

  out <- new_testflow(
    "linear_regression",
    "multiple linear regression",
    outcome_nm,
    paste(predictor_nms, collapse = ", "),
    data = df,
    descriptives = descriptives_numeric(df, c(outcome_nm, predictor_nms)),
    assumptions = assumption_checks(check_independence_note(), normality, homoscedasticity, multicollinearity),
    recommended = list(test = "Linear regression", rationale = "Continuous outcome modeled as a linear function of the predictors."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = list(coefficients = coefs),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "linear_regression"
  )
  out$interpretation <- make_report(out, alpha)
  out
}

#' Logistic regression workflow
#'
#' @param formula A formula such as `outcome ~ x1 + x2 + x3`, with a binary
#'   (two-level factor or 0/1 numeric) outcome.
#' @param data A data frame.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_logistic_regression`. The
#' object is a list containing the cleaned data, descriptive statistics,
#' multicollinearity/influential-point assumption checks, the likelihood-ratio
#' test as the primary result with null hypothesis, the per-term coefficient
#' tables on the log-odds and odds-ratio scale, McFadden's pseudo R squared as
#' the effect size, optional `ggplot`, original call, and report text.
#' @details
#' With outcome \eqn{y_i \in \{0,1\}} and predictors
#' \eqn{x_1, \ldots, x_p}:
#' \deqn{\log\frac{P(y_i=1)}{1-P(y_i=1)} = \beta_0 + \beta_1 x_{1i} + \cdots + \beta_p x_{pi}}
#'
#' Coefficients are reported on the log-odds scale and, exponentiated, as
#' odds ratios \eqn{OR_j = e^{\beta_j}}.
#'
#' The overall model test is the likelihood-ratio test against the
#' intercept-only model, using each model's residual deviance
#' \eqn{D = -2\log L}:
#' \deqn{LR = D_{null} - D_{model} \sim \chi^2_{p}}
#'
#' McFadden's pseudo R squared:
#' \deqn{R^2_{McFadden} = 1 - \frac{\log L_{model}}{\log L_{null}}}
#' @references
#' Hosmer DW, Lemeshow S, Sturdivant RX. *Applied Logistic Regression* (3rd
#' ed.). Wiley; 2013.
#' @export
test_logistic_regression <- function(formula, data, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  formula_obj <- stats::as.formula(deparse(formula))
  outcome_nm <- all.vars(formula_obj[[2]])[1]
  predictor_nms <- all.vars(formula_obj[[3]])
  df <- drop_missing(data, c(outcome_nm, predictor_nms), na.rm = na.rm)

  outcome_vals <- unique(stats::na.omit(df[[outcome_nm]]))
  if (length(outcome_vals) != 2) {
    stop("`", outcome_nm, "` must have exactly two non-missing levels for logistic regression.", call. = FALSE)
  }
  if (!is.factor(df[[outcome_nm]])) df[[outcome_nm]] <- as.factor(df[[outcome_nm]])

  fit <- stats::glm(formula_obj, data = df, family = stats::binomial())
  null_fit <- stats::update(fit, . ~ 1)

  coefs_log_odds <- broom::tidy(fit, conf.int = TRUE)
  coefs_or <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  names(coefs_or)[names(coefs_or) %in% c("estimate", "conf.low", "conf.high")] <- c("odds_ratio", "or.conf.low", "or.conf.high")

  multicollinearity <- if (length(predictor_nms) < 2) {
    assumption_check("Multicollinearity", "not applicable", "Only one predictor; multicollinearity does not apply.")
  } else if (requireNamespace("car", quietly = TRUE)) {
    vifs <- car::vif(fit)
    max_vif <- max(vifs)
    assumption_check(
      "Multicollinearity",
      ifelse(max_vif < 5, "acceptable", "not acceptable"),
      ifelse(max_vif < 5, "Variance inflation factors are all below 5.", "At least one predictor has a variance inflation factor of 5 or higher, indicating multicollinearity."),
      method = "Variance inflation factor",
      statistic = max_vif,
      details = paste0("Max VIF (", names(vifs)[which.max(vifs)], ") = ", format_stat(max_vif))
    )
  } else {
    assumption_check("Multicollinearity", "not checked", "Install the 'car' package to compute variance inflation factors.")
  }

  cooks_d <- stats::cooks.distance(fit)
  cooks_threshold <- 4 / length(cooks_d)
  n_influential <- sum(cooks_d > cooks_threshold)
  influential <- assumption_check(
    "Influential observations",
    ifelse(n_influential == 0, "acceptable", "warning"),
    ifelse(n_influential == 0, "No observations exceed the Cook's distance threshold.", paste0(n_influential, " observation(s) exceed the Cook's distance threshold (4/n); consider inspecting them.")),
    method = "Cook's distance",
    statistic = max(cooks_d),
    details = paste0("Threshold = ", format_stat(cooks_threshold, 4))
  )

  lr_stat <- fit$null.deviance - fit$deviance
  lr_df <- fit$df.null - fit$df.residual
  lr_p <- stats::pchisq(lr_stat, lr_df, lower.tail = FALSE)
  effect <- mcfadden_r2_glm(fit, null_fit)
  primary <- tibble::tibble(
    method = "Logistic regression (likelihood-ratio test)",
    statistic = lr_stat,
    parameter = lr_df,
    p.value = lr_p
  )
  h0 <- paste0("H0: none of ", paste(predictor_nms, collapse = ", "), " are associated with ", outcome_nm, " (all log-odds slopes equal 0).")

  plt <- if (plot) {
    or_df <- coefs_or[coefs_or$term != "(Intercept)", , drop = FALSE]
    ggplot2::ggplot(or_df, ggplot2::aes(x = .data$odds_ratio, y = stats::reorder(.data$term, .data$odds_ratio))) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "#8A8F98") +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = .data$or.conf.low, xmax = .data$or.conf.high), color = "#4C78A8") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Logistic regression workflow", subtitle = "Odds ratios with 95% confidence intervals", x = "Odds ratio (log scale)", y = NULL) +
      ggplot2::theme_minimal()
  } else NULL

  out <- new_testflow(
    "logistic_regression",
    "logistic regression",
    outcome_nm,
    paste(predictor_nms, collapse = ", "),
    data = df,
    descriptives = descriptives_categorical(df, outcome_nm),
    assumptions = assumption_checks(check_independence_note(), multicollinearity, influential),
    recommended = list(test = "Logistic regression", rationale = "Binary outcome modeled as a logistic function of the predictors."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = list(coefficients = coefs_log_odds, odds_ratios = coefs_or),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "logistic_regression"
  )
  out$interpretation <- make_report(out, alpha)
  out
}
