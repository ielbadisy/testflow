#' Diagnostic test accuracy workflow
#'
#' @param data A data frame.
#' @param test Binary test-result column (two levels; the alphabetically
#'   second level is treated as "positive").
#' @param reference Binary gold-standard column (two levels; the
#'   alphabetically second level is treated as "positive").
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_diagnostic`. The object
#' is a list containing the cleaned data, categorical descriptives, the
#' accuracy-vs-no-information-rate test as the primary result with null
#' hypothesis, a table of sensitivity/specificity/predictive values/
#' likelihood ratios with confidence intervals, the confusion matrix,
#' accuracy as the effect size, optional `ggplot`, original call, and report
#' text.
#' @details
#' With true positives \eqn{TP}, false positives \eqn{FP}, false negatives
#' \eqn{FN}, and true negatives \eqn{TN}:
#' \deqn{Sensitivity = \frac{TP}{TP+FN}, \qquad Specificity = \frac{TN}{TN+FP}}
#' \deqn{PPV = \frac{TP}{TP+FP}, \qquad NPV = \frac{TN}{TN+FN}}
#' \deqn{LR^+ = \frac{Sensitivity}{1-Specificity}, \qquad
#' LR^- = \frac{1-Sensitivity}{Specificity}}
#'
#' Sensitivity, specificity, PPV, NPV, and accuracy each get an exact
#' (Clopper-Pearson) confidence interval. Likelihood ratios get the
#' closed-form log-scale interval of Simel, Samsa & Matchar (1991):
#' \deqn{SE[\log LR^+] = \sqrt{\frac1{TP}-\frac1{TP+FN}+\frac1{FP}-\frac1{FP+TN}}}
#'
#' The primary test compares overall accuracy to the no-information rate
#' (the larger of the two reference-class proportions) via an exact binomial
#' test, following the same convention as `caret::confusionMatrix()`.
#' @references
#' Simel DL, Samsa GP, Matchar DB. Likelihood ratios with confidence:
#' sample size estimation for diagnostic test studies. *Journal of Clinical
#' Epidemiology*. 1991;44(8):763-770.
#'
#' Altman DG, Bland JM. Diagnostic tests 1: sensitivity and specificity.
#' *BMJ*. 1994;308(6943):1552.
#' @export
test_diagnostic <- function(data, test, reference, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  test_nm <- rlang::as_name(rlang::ensym(test))
  reference_nm <- rlang::as_name(rlang::ensym(reference))
  df <- drop_missing(data, c(test_nm, reference_nm), na.rm = na.rm)

  test_levels <- assert_two_groups(df, test_nm)
  reference_levels <- assert_two_groups(df, reference_nm)
  df[[test_nm]] <- factor(df[[test_nm]], levels = test_levels)
  df[[reference_nm]] <- factor(df[[reference_nm]], levels = reference_levels)

  tab <- table(test = df[[test_nm]], reference = df[[reference_nm]])
  tn <- unname(tab[1, 1]); fp <- unname(tab[2, 1]); fn <- unname(tab[1, 2]); tp <- unname(tab[2, 2])
  n <- tp + fp + fn + tn

  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  accuracy <- (tp + tn) / n
  lr_pos <- sensitivity / (1 - specificity)
  lr_neg <- (1 - sensitivity) / specificity

  exact_ci <- function(x, m) {
    if (m == 0) return(c(NA_real_, NA_real_))
    as.numeric(stats::binom.test(x, m, conf.level = 1 - alpha)$conf.int)
  }
  sens_ci <- exact_ci(tp, tp + fn)
  spec_ci <- exact_ci(tn, tn + fp)
  ppv_ci <- exact_ci(tp, tp + fp)
  npv_ci <- exact_ci(tn, tn + fn)
  acc_ci <- exact_ci(tp + tn, n)

  z <- stats::qnorm(1 - alpha / 2)
  se_log_lrp <- sqrt(1 / tp - 1 / (tp + fn) + 1 / fp - 1 / (fp + tn))
  lrp_ci <- exp(log(lr_pos) + c(-1, 1) * z * se_log_lrp)
  se_log_lrn <- sqrt(1 / fn - 1 / (tp + fn) + 1 / tn - 1 / (fp + tn))
  lrn_ci <- exp(log(lr_neg) + c(-1, 1) * z * se_log_lrn)

  diagnostic_table <- tibble::tibble(
    metric = c("Sensitivity", "Specificity", "Positive predictive value", "Negative predictive value", "Accuracy", "Positive likelihood ratio", "Negative likelihood ratio"),
    estimate = c(sensitivity, specificity, ppv, npv, accuracy, lr_pos, lr_neg),
    conf.low = c(sens_ci[1], spec_ci[1], ppv_ci[1], npv_ci[1], acc_ci[1], lrp_ci[1], lrn_ci[1]),
    conf.high = c(sens_ci[2], spec_ci[2], ppv_ci[2], npv_ci[2], acc_ci[2], lrp_ci[2], lrn_ci[2])
  )

  nir <- max((tp + fn) / n, (fp + tn) / n)
  acc_test <- stats::binom.test(tp + tn, n, p = nir, alternative = "greater")
  primary <- safe_tidy_htest(acc_test, "Accuracy vs. no-information rate")
  h0 <- paste0("H0: accuracy equals the no-information rate (", format_stat(nir), ").")
  effect <- tibble::tibble(name = "Accuracy", estimate = accuracy, magnitude = NA_character_)

  plt <- if (plot) {
    pd <- diagnostic_table[diagnostic_table$metric %in% c("Sensitivity", "Specificity", "Positive predictive value", "Negative predictive value", "Accuracy"), ]
    ggplot2::ggplot(pd, ggplot2::aes(x = .data$metric, y = .data$estimate)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), color = "#4C78A8") +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(formatC(100 * x, format = "f", digits = 0), "%")) +
      ggplot2::labs(title = "Diagnostic accuracy workflow", subtitle = paste0("Accuracy = ", format_stat(accuracy * 100), "%, vs. no-information rate ", format_stat(nir * 100), "%"), x = NULL, y = NULL) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal()
  } else NULL

  out <- new_testflow(
    "diagnostic",
    "diagnostic test accuracy",
    test_nm,
    reference_nm,
    data = df,
    descriptives = descriptives_categorical(df, c(test_nm, reference_nm)),
    assumptions = assumption_checks(check_independence_note("Test and reference results are assumed to be measured independently.")),
    recommended = list(test = "Diagnostic accuracy (2x2 table)", rationale = "Test result evaluated against a binary gold-standard reference."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = list(diagnostic_table = diagnostic_table, confusion_matrix = tab),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "diagnostic"
  )
  out$interpretation <- make_report(out, alpha)
  out
}

#' Receiver operating characteristic (ROC) curve workflow
#'
#' @param data A data frame.
#' @param predictor Numeric predictor/biomarker column; higher values are
#'   assumed to indicate the positive class.
#' @param outcome Binary outcome column (two levels; the alphabetically
#'   second level is treated as the positive class).
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_roc`. The object is a
#' list containing the cleaned data, descriptive statistics by outcome
#' class, the test of AUC against 0.5 as the primary result with null
#' hypothesis, the full ROC curve and the Youden's-J-optimal threshold, AUC
#' as the effect size, optional `ggplot`, original call, and report text.
#' @details
#' The area under the ROC curve equals the probability that a randomly
#' chosen positive-class observation has a higher predictor value than a
#' randomly chosen negative-class observation, and is computed from the
#' Mann-Whitney U statistic:
#' \deqn{AUC = \frac{U}{n_1 n_0}}
#'
#' The standard error uses the closed-form Hanley & McNeil (1982) formula
#' (an asymptotically equivalent, but not numerically identical, alternative
#' to DeLong's method):
#' \deqn{SE(AUC) = \sqrt{\frac{AUC(1-AUC)+(n_1-1)(Q_1-AUC^2)+(n_0-1)(Q_2-AUC^2)}{n_1 n_0}}}
#' \deqn{Q_1 = \frac{AUC}{2-AUC}, \qquad Q_2 = \frac{2AUC^2}{1+AUC}}
#'
#' The Youden's J statistic identifies the threshold maximizing
#' \eqn{J = Sensitivity + Specificity - 1}.
#' @references
#' Hanley JA, McNeil BJ. The meaning and use of the area under a receiver
#' operating characteristic (ROC) curve. *Radiology*. 1982;143(1):29-36.
#'
#' Youden WJ. Index for rating diagnostic tests. *Cancer*.
#' 1950;3(1):32-35.
#'
#' Hosmer DW, Lemeshow S, Sturdivant RX. *Applied Logistic Regression* (3rd
#' ed.). Wiley; 2013. (AUC magnitude convention, p. 177.)
#' @export
test_roc <- function(data, predictor, outcome, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  predictor_nm <- rlang::as_name(rlang::ensym(predictor))
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  df <- drop_missing(data, c(predictor_nm, outcome_nm), na.rm = na.rm)

  outcome_levels <- assert_two_groups(df, outcome_nm)
  df[[outcome_nm]] <- factor(df[[outcome_nm]], levels = outcome_levels)
  pos <- df[[predictor_nm]][df[[outcome_nm]] == outcome_levels[2]]
  neg <- df[[predictor_nm]][df[[outcome_nm]] == outcome_levels[1]]
  n_pos <- length(pos); n_neg <- length(neg)

  wt <- stats::wilcox.test(pos, neg, exact = FALSE)
  auc <- unname(wt$statistic) / (n_pos * n_neg)

  q1 <- auc / (2 - auc)
  q2 <- 2 * auc^2 / (1 + auc)
  se_auc <- sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) + (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg))
  z <- stats::qnorm(1 - alpha / 2)
  auc_ci <- auc + c(-1, 1) * z * se_auc
  z_stat <- (auc - 0.5) / se_auc
  p_value <- 2 * stats::pnorm(-abs(z_stat))

  ord <- order(-df[[predictor_nm]])
  y_sorted <- df[[outcome_nm]][ord] == outcome_levels[2]
  pred_sorted <- df[[predictor_nm]][ord]
  cum_tp <- cumsum(y_sorted)
  cum_fp <- cumsum(!y_sorted)
  last_idx <- !duplicated(pred_sorted, fromLast = TRUE)
  roc_curve <- tibble::tibble(
    threshold = c(Inf, pred_sorted[last_idx]),
    sensitivity = c(0, cum_tp[last_idx] / n_pos),
    specificity = c(1, 1 - cum_fp[last_idx] / n_neg)
  )
  youden <- roc_curve$sensitivity + roc_curve$specificity - 1
  best <- which.max(youden)
  optimal <- roc_curve[best, ]

  effect <- tibble::tibble(name = "AUC", estimate = auc, magnitude = magnitude_auc(auc))
  primary <- tibble::tibble(
    method = "AUC vs. 0.5 (no discrimination)",
    statistic = z_stat,
    p.value = p_value,
    conf.low = auc_ci[1],
    conf.high = auc_ci[2]
  )
  h0 <- paste0("H0: AUC = 0.5 (", predictor_nm, " does not discriminate between levels of ", outcome_nm, ").")

  plt <- if (plot) {
    ggplot2::ggplot(roc_curve, ggplot2::aes(x = 1 - .data$specificity, y = .data$sensitivity)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#8A8F98") +
      ggplot2::geom_step(color = "#4C78A8", linewidth = 0.8) +
      ggplot2::geom_point(data = optimal, ggplot2::aes(x = 1 - .data$specificity, y = .data$sensitivity), color = "#F58518", size = 2.5) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = "ROC curve workflow",
        subtitle = paste0("AUC = ", format_stat(auc), " (", effect$magnitude[1], ")"),
        x = "1 - Specificity", y = "Sensitivity"
      ) +
      ggplot2::theme_minimal()
  } else NULL

  out <- new_testflow(
    "roc",
    "ROC curve",
    predictor_nm,
    outcome_nm,
    data = df,
    descriptives = descriptives_numeric(df, predictor_nm, outcome_nm),
    assumptions = assumption_checks(check_independence_note("Higher predictor values are assumed to indicate the positive class.")),
    recommended = list(test = "ROC / AUC analysis", rationale = "Continuous predictor evaluated against a binary outcome."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = list(
      roc_curve = roc_curve,
      optimal_threshold = tibble::tibble(threshold = optimal$threshold, sensitivity = optimal$sensitivity, specificity = optimal$specificity, youden_j = youden[best])
    ),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "roc"
  )
  out$interpretation <- make_report(out, alpha)
  out
}
