#' Inter-rater agreement workflow (Cohen's kappa)
#'
#' @param data A data frame.
#' @param rater1 First rater's categorical column.
#' @param rater2 Second rater's categorical column, using the same category
#'   set as `rater1`.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_agreement`. The object
#' is a list containing the cleaned data, categorical descriptives, the test
#' of kappa against 0 as the primary result with null hypothesis, the
#' agreement (confusion) table, Cohen's kappa as the effect size, optional
#' `ggplot`, original call, and report text.
#' @details
#' With observed agreement \eqn{p_o} (the proportion of subjects on which
#' both raters agree) and chance-expected agreement \eqn{p_e} (from the
#' marginal category proportions):
#' \deqn{\kappa = \frac{p_o - p_e}{1 - p_e}}
#'
#' Two different standard errors are used for two different purposes, since
#' the sampling variance of \eqn{\hat\kappa} depends on \eqn{\kappa} itself:
#'
#' - The confidence interval uses the large-sample SE of Fleiss, Cohen &
#'   Everitt (1969), evaluated *at* \eqn{\hat\kappa}, which accounts for the
#'   full marginal covariance structure (not the simpler
#'   \eqn{\sqrt{p_o(1-p_o)}/[(1-p_e)\sqrt n]} approximation sometimes seen,
#'   which understates the variance when categories are unevenly
#'   distributed).
#' - The z-test of \eqn{H_0:\kappa=0} uses a *different*, smaller-magnitude
#'   SE (Fleiss, 1981), evaluated under the null \eqn{\kappa=0} rather than
#'   at \eqn{\hat\kappa}: \eqn{SE_0=\sqrt{p_e+p_e^2-\sum_i\pi_{i\cdot}\pi_{\cdot
#'   i}(\pi_{i\cdot}+\pi_{\cdot i})}/[(1-p_e)\sqrt n]}. Using the
#'   confidence-interval SE for this test instead (a mistake this function
#'   previously made) understates the null SE whenever agreement is
#'   materially above chance - the typical case - inflating the z-statistic
#'   and anti-conservatively shrinking the p-value; verified in package
#'   tests to match `irr::kappa2()`'s reported z exactly.
#' @references
#' Cohen J. A coefficient of agreement for nominal scales. *Educational and
#' Psychological Measurement*. 1960;20(1):37-46.
#'
#' Fleiss JL, Cohen J, Everitt BS. Large sample standard errors of kappa and
#' weighted kappa. *Psychological Bulletin*. 1969;72(5):323-327.
#'
#' Fleiss JL. *Statistical Methods for Rates and Proportions* (2nd ed.).
#' Wiley; 1981. (Null-hypothesis SE for testing kappa = 0.)
#'
#' Landis JR, Koch GG. The measurement of observer agreement for
#' categorical data. *Biometrics*. 1977;33(1):159-174. (Magnitude
#' convention.)
#' @export
test_agreement <- function(data, rater1, rater2, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  rater1_nm <- rlang::as_name(rlang::ensym(rater1))
  rater2_nm <- rlang::as_name(rlang::ensym(rater2))
  df <- drop_missing(data, c(rater1_nm, rater2_nm), na.rm = na.rm)

  category_levels <- sort(unique(c(as.character(df[[rater1_nm]]), as.character(df[[rater2_nm]]))))
  df[[rater1_nm]] <- factor(df[[rater1_nm]], levels = category_levels)
  df[[rater2_nm]] <- factor(df[[rater2_nm]], levels = category_levels)
  tab <- table(df[[rater1_nm]], df[[rater2_nm]])

  n <- sum(tab)
  p <- tab / n
  pi_row <- rowSums(p)
  pi_col <- colSums(p)
  po <- sum(diag(p))
  pe <- sum(pi_row * pi_col)
  kappa <- (po - pe) / (1 - pe)

  k <- nrow(p)
  term1 <- sum(vapply(seq_len(k), function(i) p[i, i] * (1 - (pi_row[i] + pi_col[i]) * (1 - kappa))^2, numeric(1)))
  term2 <- (1 - kappa)^2 * sum(vapply(seq_len(k), function(i) {
    sum(vapply(seq_len(k), function(j) if (i == j) 0 else p[i, j] * (pi_col[i] + pi_row[j])^2, numeric(1)))
  }, numeric(1)))
  term3 <- (kappa - pe * (1 - kappa))^2
  se_kappa <- sqrt((term1 + term2 - term3) / (n * (1 - pe)^2))

  # Null-hypothesis SE (Fleiss, 1981), evaluated at kappa = 0 rather than at
  # the observed kappa_hat, for the z-test of H0: kappa = 0. This is
  # deliberately a *different* SE from se_kappa above: se_kappa (Fleiss,
  # Cohen & Everitt 1969) is the asymptotic SE of kappa_hat around its own
  # true value and is the correct SE for a confidence interval, but reusing
  # it as the denominator of a null test is a known error (it understates
  # se0_kappa whenever agreement is not near chance, which happens to be the
  # typical case, inflating the z-statistic and anti-conservatively
  # shrinking the p-value). Matches `irr::kappa2()`'s reported z exactly.
  se0_kappa <- sqrt(pe + pe^2 - sum(pi_row * pi_col * (pi_row + pi_col))) / ((1 - pe) * sqrt(n))

  z <- stats::qnorm(1 - alpha / 2)
  kappa_ci <- kappa + c(-1, 1) * z * se_kappa
  z_stat <- kappa / se0_kappa
  p_value <- 2 * stats::pnorm(-abs(z_stat))

  effect <- tibble::tibble(name = "Cohen's kappa", estimate = kappa, magnitude = magnitude_kappa(kappa))
  primary <- tibble::tibble(
    method = "Cohen's kappa (vs. 0)",
    statistic = z_stat,
    p.value = p_value,
    conf.low = kappa_ci[1],
    conf.high = kappa_ci[2]
  )
  h0 <- "H0: agreement beyond chance is zero (kappa = 0)."

  plt <- if (plot) {
    hm <- as.data.frame(tab)
    names(hm) <- c(rater1_nm, rater2_nm, "n")
    ggplot2::ggplot(hm, ggplot2::aes(x = .data[[rater1_nm]], y = .data[[rater2_nm]], fill = .data$n)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = .data$n), color = "white") +
      ggplot2::scale_fill_gradient(low = "#A9BEDB", high = "#1F2933") +
      ggplot2::coord_equal() +
      ggplot2::labs(title = "Inter-rater agreement workflow", subtitle = paste0("Cohen's kappa = ", format_stat(kappa), " (", effect$magnitude[1], ")"), x = rater1_nm, y = rater2_nm, fill = "n") +
      ggplot2::theme_minimal()
  } else NULL

  out <- new_testflow(
    "agreement",
    "inter-rater agreement",
    rater1_nm,
    rater2_nm,
    data = df,
    descriptives = descriptives_categorical(df, c(rater1_nm, rater2_nm)),
    assumptions = assumption_checks(check_independence_note("Raters are assumed to classify subjects independently.")),
    recommended = list(test = "Cohen's kappa", rationale = "Two raters classifying the same subjects into the same categories."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = list(agreement_table = tab),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "agreement"
  )
  out$interpretation <- make_report(out, alpha)
  out
}

#' Intraclass correlation coefficient workflow
#'
#' @param data A data frame in wide format (one column per rater/measurement,
#'   one row per subject).
#' @param measures Rater/measurement columns selected with tidyselect syntax.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_icc`. The object is a
#' list containing the cleaned long-format data, descriptive statistics, the
#' F-test for ICC(2,1) as the primary result with null hypothesis, a table
#' of the one-way (ICC1), two-way random/absolute-agreement (ICC2), and
#' two-way fixed/consistency (ICC3) estimates with confidence intervals,
#' ICC(2,1) as the effect size, optional `ggplot`, original call, and report
#' text.
#' @details
#' Following Shrout & Fleiss (1979), variance components come from one-way
#' (\eqn{BMS}, \eqn{WMS}) and two-way (\eqn{BMS}, \eqn{JMS}, \eqn{EMS})
#' ANOVA decompositions of the \eqn{n} subjects by \eqn{k} raters:
#' \deqn{ICC(1,1) = \frac{BMS-WMS}{BMS+(k-1)WMS}}
#' \deqn{ICC(2,1) = \frac{BMS-EMS}{BMS+(k-1)EMS+\frac{k}{n}(JMS-EMS)}
#' \quad\text{(two-way random, absolute agreement)}}
#' \deqn{ICC(3,1) = \frac{BMS-EMS}{BMS+(k-1)EMS}
#' \quad\text{(two-way fixed, consistency)}}
#'
#' ICC(2,1) (absolute agreement, allowing raters to differ systematically)
#' is reported as the primary/effect-size estimate, following the
#' single-measurement, absolute-agreement, two-way random-effects
#' recommendation of Koo & Li (2016) for typical reliability studies.
#' Confidence intervals use the F-distribution-based formulas of McGraw &
#' Wong (1996); ICC(2,1)'s interval uses a Satterthwaite-approximated
#' denominator degrees of freedom.
#' @references
#' Shrout PE, Fleiss JL. Intraclass correlations: uses in assessing rater
#' reliability. *Psychological Bulletin*. 1979;86(2):420-428.
#'
#' McGraw KO, Wong SP. Forming inferences about some intraclass correlation
#' coefficients. *Psychological Methods*. 1996;1(1):30-46.
#'
#' Koo TK, Li MY. A guideline of selecting and reporting intraclass
#' correlation coefficients for reliability research. *Journal of
#' Chiropractic Medicine*. 2016;15(2):155-163.
#' @export
test_icc <- function(data, measures, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  measure_nms <- tidyselect_names(data, {{ measures }})
  df <- drop_missing(data, measure_nms, na.rm = na.rm)
  df$.subject <- factor(seq_len(nrow(df)))
  long <- tidyr::pivot_longer(df, dplyr::all_of(measure_nms), names_to = ".rater", values_to = ".value")
  long$.rater <- factor(long$.rater)

  n <- nrow(df)
  k <- length(measure_nms)
  z <- stats::qnorm(1 - alpha / 2)

  fit1 <- stats::aov(.value ~ .subject, data = long)
  tab1 <- stats::anova(fit1)
  bms <- tab1[".subject", "Mean Sq"]
  wms <- tab1["Residuals", "Mean Sq"]
  icc1 <- (bms - wms) / (bms + (k - 1) * wms)
  f1 <- bms / wms
  fl1 <- f1 / stats::qf(1 - alpha / 2, n - 1, n * (k - 1))
  fu1 <- f1 * stats::qf(1 - alpha / 2, n * (k - 1), n - 1)
  icc1_ci <- c((fl1 - 1) / (fl1 + k - 1), (fu1 - 1) / (fu1 + k - 1))
  icc1_p <- stats::pf(f1, n - 1, n * (k - 1), lower.tail = FALSE)

  fit2 <- stats::aov(.value ~ .subject + .rater, data = long)
  tab2 <- stats::anova(fit2)
  bms2 <- tab2[".subject", "Mean Sq"]
  jms <- tab2[".rater", "Mean Sq"]
  ems <- tab2["Residuals", "Mean Sq"]

  icc3 <- (bms2 - ems) / (bms2 + (k - 1) * ems)
  f3 <- bms2 / ems
  df2_3 <- (n - 1) * (k - 1)
  fl3 <- f3 / stats::qf(1 - alpha / 2, n - 1, df2_3)
  fu3 <- f3 * stats::qf(1 - alpha / 2, df2_3, n - 1)
  icc3_ci <- c((fl3 - 1) / (fl3 + k - 1), (fu3 - 1) / (fu3 + k - 1))
  icc3_p <- stats::pf(f3, n - 1, df2_3, lower.tail = FALSE)

  icc2 <- (bms2 - ems) / (bms2 + (k - 1) * ems + (k / n) * (jms - ems))
  # The F-test (H0: ICC = 0) and the confidence interval use different
  # Satterthwaite coefficients: the test evaluates a/b at the null value
  # r0 = 0 (which collapses to a = 0, b = 1, i.e. F2 = BMS/EMS under the
  # null), while the CI evaluates a/b at the *estimated* icc2 - these are
  # not the same computation and must not share one a/b pair.
  f2 <- bms2 / ems
  df2_test <- (n - 1) * (k - 1)
  icc2_p <- stats::pf(f2, n - 1, df2_test, lower.tail = FALSE)

  a_coef <- (k * icc2) / (n * (1 - icc2))
  b_coef <- 1 + (k * icc2 * (n - 1)) / (n * (1 - icc2))
  v_df <- (a_coef * jms + b_coef * ems)^2 / ((a_coef * jms)^2 / (k - 1) + (b_coef * ems)^2 / ((n - 1) * (k - 1)))
  fl2 <- stats::qf(1 - alpha / 2, n - 1, v_df)
  fu2 <- stats::qf(1 - alpha / 2, v_df, n - 1)
  icc2_ci <- c(
    (n * (bms2 - fl2 * ems)) / (fl2 * (k * jms + (k * n - k - n) * ems) + n * bms2),
    (n * (fu2 * bms2 - ems)) / (k * jms + (k * n - k - n) * ems + n * fu2 * bms2)
  )

  icc_table <- tibble::tibble(
    method = c("ICC(1,1) one-way random", "ICC(2,1) two-way random, absolute agreement", "ICC(3,1) two-way fixed, consistency"),
    estimate = c(icc1, icc2, icc3),
    conf.low = c(icc1_ci[1], icc2_ci[1], icc3_ci[1]),
    conf.high = c(icc1_ci[2], icc2_ci[2], icc3_ci[2]),
    statistic = c(f1, f2, f3),
    p.value = c(icc1_p, icc2_p, icc3_p)
  )

  effect <- tibble::tibble(name = "ICC(2,1)", estimate = icc2, magnitude = magnitude_icc(icc2))
  primary <- tibble::tibble(method = "ICC(2,1) F-test (vs. 0)", statistic = f2, parameter = n - 1, p.value = icc2_p, conf.low = icc2_ci[1], conf.high = icc2_ci[2])
  h0 <- "H0: ICC(2,1) = 0 (no reliability beyond chance)."

  plt <- if (plot) {
    ggplot2::ggplot(icc_table, ggplot2::aes(x = .data$method, y = .data$estimate)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high), color = "#4C78A8") +
      ggplot2::scale_y_continuous(limits = c(min(0, icc_table$conf.low, na.rm = TRUE), 1)) +
      ggplot2::labs(title = "Intraclass correlation workflow", subtitle = paste0("ICC(2,1) = ", format_stat(icc2), " (", effect$magnitude[1], ")"), x = NULL, y = "ICC") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal()
  } else NULL

  out <- new_testflow(
    "icc",
    "intraclass correlation",
    paste(measure_nms, collapse = ", "),
    data = long,
    descriptives = descriptives_numeric(long, ".value", ".rater"),
    assumptions = assumption_checks(check_independence_note("Subjects are assumed to be measured independently of one another.")),
    recommended = list(test = "Intraclass correlation coefficient", rationale = "Multiple raters or repeated measurements of the same subjects."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = list(icc_table = icc_table),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "icc"
  )
  out$interpretation <- make_report(out, alpha)
  out
}
