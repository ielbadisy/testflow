#' Kaplan-Meier and log-rank survival workflow
#'
#' @param formula A formula such as `Surv(time, status) ~ group`, with a
#'   two-level grouping factor.
#' @param data A data frame.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param risktable Logical; if `TRUE` (default) the numbers at risk per
#'   group are printed under the curves at the x-axis breaks.
#' @param times Optional x-axis break times (and risk-table times). Default
#'   `NULL` uses pretty breaks over the follow-up.
#' @return A `testflow` object with class `testflow_survival`. The object is
#' a list containing the cleaned data, per-group descriptive survival
#' statistics, the log-rank test as the primary result with null hypothesis,
#' the Kaplan-Meier curve data, a companion univariate hazard ratio as the
#' effect size, optional `ggplot`, original call, and report text. The
#' log-rank p-value is in `$primary_test$p.value` (there is no top-level
#' `$p_value`). The hazard ratio, with its 95% Wald CI (`conf.low`,
#' `conf.high`), is in `$effect_size`; its name states the comparison, for
#' example "Hazard ratio (F vs M)", where the second factor level is compared
#' to the first (reference) level. Relevel the factor to change the reference.
#' The Cox model table is in `$alternative_tests$cox`.
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
testsurvival <- function(formula, data, alpha = 0.05, plot = TRUE, na.rm = TRUE, risktable = TRUE, times = NULL) {
  formula_obj <- stats::as.formula(paste(deparse(formula), collapse = " "))
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
  descriptives <- .bt(
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
  lv <- levels(df[[group_nm]])
  # The Cox HR compares the second level to the first (reference) level;
  # name both so the direction is explicit when level order changes.
  effect <- .bt(
    name = paste0("Hazard ratio (", lv[2], " vs ", lv[1], ")"),
    estimate = hr,
    magnitude = NA_character_,
    conf.low = unname(hr_ci[1]),
    conf.high = unname(hr_ci[2])
  )

  # The log-rank test has no confidence interval; the HR CI lives with the
  # effect size and the Cox table, not on the log-rank row.
  primary <- .bt(
    method = "Log-rank test",
    statistic = chisq,
    parameter = lr_df,
    p.value = lr_p
  )
  h0 <- paste0("H0: the survival distributions are equal across levels of ", group_nm, ".")

  plt <- if (plot) km_plot(fit_km, lv, group_nm, time_nm, risktable, times,
                           paste0("Log-rank test, ", pstr(lr_p), "; ", effstr(effect))) else NULL

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
#' @param units Optional named numeric vector of plotting units for numeric
#'   predictors, for example `c(age = 10)` to plot the hazard ratio per 10
#'   years. Only the plot changes; tables stay per unit.
#' @param crude Logical; if `TRUE`, the plot adds the crude hazard ratio of
#'   each predictor (one univariable Cox model per predictor, same rows) in
#'   grey under the adjusted one. The crude table is stored in
#'   `$alternative_tests$crude`. Default `FALSE`.
#' @return A `testflow` object with class `testflow_cox`. The object is a
#' list containing the cleaned data, descriptive statistics, proportional-
#' hazards assumption checks (the global `cox.zph()` test, then one row per
#' term named `"Proportional hazards (<term>)"`), the overall likelihood-ratio test as the
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
testcox <- function(formula, data, alpha = 0.05, plot = TRUE, na.rm = TRUE, units = NULL, crude = FALSE) {
  formula_obj <- stats::as.formula(paste(deparse(formula), collapse = " "))
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
  zt <- zph$table[rownames(zph$table) != "GLOBAL", , drop = FALSE]
  ph_terms <- .map_dfr(rownames(zt), function(tm) assumption_check(
    paste0("Proportional hazards (", tm, ")"), ifelse(zt[tm, "p"] >= alpha, "acceptable", "not acceptable"),
    ifelse(zt[tm, "p"] >= alpha, "No evidence that this term's effect changes over time.",
           "This term's scaled Schoenfeld residuals trend with time; its hazard ratio may vary over follow-up."),
    method = "Schoenfeld residual test (Grambsch-Therneau)", statistic = unname(zt[tm, "chisq"]), p_value = unname(zt[tm, "p"]),
    details = paste0("df = ", zt[tm, "df"])))

  effect <- .bt(name = "Concordance index", estimate = unname(g$concordance), magnitude = NA_character_)
  primary <- .bt(
    method = "Cox proportional hazards regression (likelihood-ratio test)",
    statistic = unname(logtest["test"]),
    parameter = unname(logtest["df"]),
    p.value = unname(logtest["pvalue"])
  )
  h0 <- paste0("H0: none of ", paste(predictor_nms, collapse = ", "), " are associated with the hazard (all log hazard ratios equal 0).")

  cr <- if (crude) .map_dfr(attr(stats::terms(fit), "term.labels"), function(tm) broom::tidy(survival::coxph(
    stats::reformulate(tm, formula_obj[[2]]), data = df), conf.int = TRUE))
  plt <- if (plot) ratio_plot(coefs_log, "Hazard ratio, 95% CI (log scale)", predictor_nms, units, crude = cr,
                              subtitle = paste0("Likelihood-ratio test, ", pstr(unname(logtest["pvalue"])))) else NULL

  out <- new_testflow(
    "cox",
    "Cox proportional hazards regression",
    time_nm,
    paste(predictor_nms, collapse = ", "),
    data = df,
    descriptives = descriptives_numeric(df, intersect(predictor_nms, names(df)[vapply(df, is.numeric, logical(1))])),
    assumptions = assumption_checks(check_independence_note("Independent censoring is assumed."), ph_check, ph_terms),
    recommended = list(test = "Cox proportional hazards regression", rationale = "Time-to-event outcome modeled as a function of the predictors."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = c(list(coefficients = coefs_log, hazard_ratios = coefs_hr), if (crude) list(crude = cr)),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "cox"
  )
  out$interpretation <- make_report(out, alpha)
  out
}

# Kaplan-Meier curves with stepped CI bands from (0, 1) and, below the axis, the numbers at risk per group.
km_plot <- function(fit, lv, group, xlab, risktable = TRUE, times = NULL, subtitle = NULL) {
  k <- length(lv); s <- summary(fit, censored = TRUE); g <- factor(sub("^.*=", "", s$strata), lv)
  d <- data.frame(time = s$time, surv = s$surv, lo = s$lower, hi = s$upper, g = g)
  d <- rbind(data.frame(time = 0, surv = 1, lo = 1, hi = 1, g = factor(lv, lv)), d)
  d <- d[order(d$g, d$time), ]
  # Stepped band: each row also ends where the next one starts.
  nx <- stats::ave(d$time, d$g, FUN = function(t) c(t[-1], t[length(t)]))
  band <- rbind(d, transform(d, time = nx)); band <- band[order(band$g, band$time, -band$surv), ]
  br <- times %||% pretty(c(0, max(d$time)))
  br <- br[br <= max(d$time)]
  col <- tfpal(k)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$time, colour = .data$g, fill = .data$g)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), data = band, alpha = 0.12, colour = NA) +
    ggplot2::geom_step(ggplot2::aes(y = .data$surv), linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = col) + ggplot2::scale_fill_manual(values = col) +
    ggplot2::scale_x_continuous(breaks = br) +
    ggplot2::labs(subtitle = subtitle, x = xlab, y = "Survival probability", colour = group, fill = group) +
    tftheme() + ggplot2::theme(legend.position = "top")
  if (!risktable) return(p + ggplot2::scale_y_continuous(limits = c(0, 1)))
  nr <- summary(fit, times = br, extend = TRUE)
  rt <- data.frame(time = nr$time, n = nr$n.risk, g = factor(sub("^.*=", "", nr$strata), lv))
  yr <- -0.09 - 0.09 * seq_len(k); rt$y <- yr[as.integer(rt$g)]
  # Group names and the "At risk" header sit left of the panel, under the y-axis labels.
  lab <- data.frame(y = c(-0.09, yr), label = c("At risk", lv))
  p + ggplot2::geom_text(ggplot2::aes(y = .data$y, label = .data$n), data = rt, size = 3, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(x = -Inf, y = .data$y, label = .data$label), data = lab, inherit.aes = FALSE,
                       hjust = 1.1, size = 3, colour = c(TFCOL[["pt"]], col)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(min(yr) - 0.04, 1)) +
    ggplot2::coord_cartesian(clip = "off")
}
