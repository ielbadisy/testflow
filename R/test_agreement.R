#' Inter-rater agreement workflow (Cohen's kappa)
#'
#' @param data A data frame.
#' @param rater1 First rater's categorical column.
#' @param rater2 Second rater's categorical column, using the same category
#'   set as `rater1`.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param weights `"unweighted"` (default) gives Cohen's kappa. `"linear"`
#'   or `"quadratic"` give Cohen's (1968) weighted kappa for ordered
#'   categories, with agreement weights \eqn{w_{ij} = 1 - |i-j|/(k-1)} or
#'   \eqn{1 - (i-j)^2/(k-1)^2}. Category order follows the factor levels
#'   (or numeric order), so code ordinal ratings as ordered factors.
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
#'
#' Weighted kappa replaces agreement on the diagonal with weighted agreement,
#' \eqn{\kappa_w = (p_{o,w} - p_{e,w})/(1 - p_{e,w})} with
#' \eqn{p_{o,w} = \sum w_{ij} p_{ij}} and
#' \eqn{p_{e,w} = \sum w_{ij} \pi_{i\cdot}\pi_{\cdot j}}. Its CI and null
#' SEs are the weighted forms in Fleiss, Cohen & Everitt (1969). Quadratic
#' weighted kappa is close to the ICC for ordinal ratings.
#' @references
#' Cohen J. Weighted kappa: nominal scale agreement with provision for
#' scaled disagreement or partial credit. *Psychological Bulletin*.
#' 1968;70(4):213-220.
#'
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
testagreement <- function(data, rater1, rater2, alpha = 0.05, plot = TRUE, na.rm = TRUE, weights = c("unweighted", "linear", "quadratic")) {
  weights <- match.arg(weights)
  rater1_nm <- rlang::as_name(rlang::ensym(rater1))
  rater2_nm <- rlang::as_name(rlang::ensym(rater2))
  df <- drop_missing(data, c(rater1_nm, rater2_nm), na.rm = na.rm)

  category_levels <- sort(unique(c(as.character(df[[rater1_nm]]), as.character(df[[rater2_nm]]))))
  if (weights != "unweighted") {
    # Weighted kappa needs the category order: factor levels, else numeric
    # order for numeric codes, else alphabetical.
    lv <- unique(c(levels(df[[rater1_nm]]), levels(df[[rater2_nm]])))
    if (length(lv)) {
      category_levels <- c(intersect(lv, category_levels), setdiff(category_levels, lv))
    } else if (is.numeric(df[[rater1_nm]]) && is.numeric(df[[rater2_nm]])) {
      category_levels <- as.character(sort(unique(c(df[[rater1_nm]], df[[rater2_nm]]))))
    }
  }
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

  kappa_name <- "Cohen's kappa"
  if (weights != "unweighted") {
    w <- agreement_weights(k, weights)
    po_w <- sum(w * p)
    pe_w <- sum(w * outer(pi_row, pi_col))
    kappa <- (po_w - pe_w) / (1 - pe_w)
    w_row <- as.vector(w %*% pi_col)
    w_col <- as.vector(crossprod(w, pi_row))
    dev_hat <- w - outer(w_row, w_col, "+") * (1 - kappa)
    se_kappa <- sqrt(max(0, sum(p * dev_hat^2) - (kappa - pe_w * (1 - kappa))^2) / (n * (1 - pe_w)^2))
    dev_0 <- w - outer(w_row, w_col, "+")
    se0_kappa <- sqrt(max(0, sum(outer(pi_row, pi_col) * dev_0^2) - pe_w^2) / (n * (1 - pe_w)^2))
    kappa_name <- paste0("Weighted kappa (", weights, ")")
  }

  z <- stats::qnorm(1 - alpha / 2)
  # Wald interval clipped to the parameter range (kappa <= 1; floor at -1).
  kappa_ci <- pmin(pmax(kappa + c(-1, 1) * z * se_kappa, -1), 1)
  z_stat <- kappa / se0_kappa
  p_value <- 2 * stats::pnorm(-abs(z_stat))

  effect <- .bt(name = kappa_name, estimate = kappa, magnitude = magnitude_kappa(kappa))
  primary <- .bt(
    method = paste0(kappa_name, " (vs. 0)"),
    statistic = z_stat,
    p.value = p_value,
    conf.low = kappa_ci[1],
    conf.high = kappa_ci[2]
  )
  h0 <- "H0: agreement beyond chance is zero (kappa = 0)."

  plt <- if (plot) {
    hm <- as.data.frame(tab, stringsAsFactors = FALSE); names(hm) <- c("r1", "r2", "n")
    lv <- rownames(tab); wr <- wraplab(12)
    hm$r1 <- factor(hm$r1, lv); hm$r2 <- factor(hm$r2, lv); hm$diag <- as.character(hm$r1) == as.character(hm$r2)
    hm$sh <- ifelse(hm$diag, "d", ifelse(hm$n > 0, "o", "z"))
    ggplot2::ggplot(hm, ggplot2::aes(x = .data$r1, y = .data$r2, fill = .data$sh)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = .data$n, colour = .data$diag), size = 3.6) +
      ggplot2::scale_fill_manual(values = c(d = TFCOL[["a1"]], o = "grey80", z = "grey95"), guide = "none") +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
      ggplot2::scale_y_discrete(labels = wr) +
      ggplot2::coord_equal() +
      ggplot2::labs(subtitle = paste0(kappa_name, " = ", format_stat(kappa), " (", 100 * (1 - alpha), "% CI ",
                                      format_stat(kappa_ci[1]), " to ", format_stat(kappa_ci[2]), ")"),
                    caption = paste0("Blue diagonal: agreement, ", sum(diag(tab)), " of ", sum(tab), " (",
                                     round(100 * sum(diag(tab)) / sum(tab)), "%)."),
                    x = rater1_nm, y = rater2_nm) +
      tftheme() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                 axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1))
  } else NULL

  out <- new_testflow(
    "agreement",
    "inter-rater agreement",
    rater1_nm,
    rater2_nm,
    data = df,
    descriptives = descriptives_categorical(df, c(rater1_nm, rater2_nm)),
    assumptions = assumption_checks(check_independence_note("Raters are assumed to classify subjects independently.")),
    recommended = list(test = kappa_name, rationale = if (weights == "unweighted") "Two raters classifying the same subjects into the same categories." else paste0("Two raters classifying the same subjects on an ordered scale; ", weights, " weights give partial credit to near-miss disagreements.")),
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

agreement_weights <- function(k, type) {
  if (k < 2) return(matrix(1, 1, 1))
  d <- abs(outer(seq_len(k), seq_len(k), "-")) / (k - 1)
  if (identical(type, "quadratic")) 1 - d^2 else 1 - d
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
testicc <- function(data, measures, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  measure_nms <- tidyselect_names(data, {{ measures }})
  df <- drop_missing(data, measure_nms, na.rm = na.rm)
  df$.subject <- factor(seq_len(nrow(df)))
  long <- .pivot_longer_simple(df, measure_nms, ".rater", ".value")
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

  icc_table <- .bt(
    method = c("ICC(1,1) one-way random", "ICC(2,1) two-way random, absolute agreement", "ICC(3,1) two-way fixed, consistency"),
    estimate = c(icc1, icc2, icc3),
    conf.low = c(icc1_ci[1], icc2_ci[1], icc3_ci[1]),
    conf.high = c(icc1_ci[2], icc2_ci[2], icc3_ci[2]),
    statistic = c(f1, f2, f3),
    p.value = c(icc1_p, icc2_p, icc3_p)
  )

  effect <- .bt(name = "ICC(2,1)", estimate = icc2, magnitude = magnitude_icc(icc2))
  primary <- .bt(method = "ICC(2,1) F-test (vs. 0)", statistic = f2, parameter = n - 1, p.value = icc2_p, conf.low = icc2_ci[1], conf.high = icc2_ci[2])
  h0 <- "H0: ICC(2,1) = 0 (no reliability beyond chance)."

  plt <- if (plot) {
    lb <- c("ICC(1,1)\none-way random", "ICC(2,1)\nabsolute agreement", "ICC(3,1)\nconsistency")
    pd <- transform(icc_table, lab = factor(lb, rev(lb)), key = c(FALSE, TRUE, FALSE))
    ggplot2::ggplot(pd, ggplot2::aes(x = .data$estimate, y = .data$lab, colour = .data$key)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = .data$conf.low, xmax = .data$conf.high)) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = TFCOL[["a1"]], `FALSE` = TFCOL[["ink"]]), guide = "none") +
      ggplot2::scale_x_continuous(limits = c(min(0, icc_table$conf.low, na.rm = TRUE), 1)) +
      ggplot2::labs(subtitle = paste0("ICC(2,1) = ", format_stat(icc2), " (", 100 * (1 - alpha), "% CI ",
                                      format_stat(icc2_ci[1]), " to ", format_stat(icc2_ci[2]), "), ", pstr(icc2_p)),
                    caption = paste0(n, " subjects, ", k, " raters."), x = "ICC", y = NULL) +
      tftheme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
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
