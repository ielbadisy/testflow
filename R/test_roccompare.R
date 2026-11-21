#' Compare two correlated ROC curves (paired DeLong test)
#'
#' Compares the AUCs of two markers (or two risk scores) measured on the
#' same subjects, using the DeLong, DeLong and Clarke-Pearson (1988)
#' nonparametric covariance of the two correlated AUCs.
#'
#' @param data A data frame with one row per subject.
#' @param predictor1,predictor2 Numeric marker columns; higher values are
#'   assumed to indicate the positive class for both.
#' @param outcome Binary outcome column (two levels; the second factor level,
#'   alphabetical unless a factor, is the positive class).
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object with both ROC curves.
#' @param na.rm Logical; remove rows with a missing value in any of the three
#'   columns (complete cases, so both AUCs use the same subjects).
#' @return A `testflow` object with class `testflow_roccompare`. The primary
#' result is the z-test of AUC1 - AUC2 = 0 with the CI of the difference.
#' `$alternative_tests$auc_table` holds each AUC with its DeLong CI and the
#' correlation between the two AUC estimates.
#' @details
#' For marker \eqn{r}, the DeLong placement values are
#' \eqn{V_{10}^r(X_i) = \frac{1}{n_0}\sum_j \psi(X_i^r, Y_j^r)} for each
#' positive \eqn{i} and \eqn{V_{01}^r(Y_j) = \frac{1}{n_1}\sum_i \psi(X_i^r, Y_j^r)}
#' for each negative \eqn{j}, with \eqn{\psi = 1} if \eqn{X > Y}, 1/2 for
#' ties, 0 otherwise. With \eqn{S_{10}} and \eqn{S_{01}} their 2 x 2 sample
#' covariance matrices,
#' \deqn{Var(\hat A_1 - \hat A_2) = \frac{c^\top S_{10} c}{n_1} + \frac{c^\top S_{01} c}{n_0}, \quad c = (1, -1)^\top}
#' and \eqn{z = (\hat A_1 - \hat A_2)/\sqrt{Var}}. This is the default paired
#' test of `pROC::roc.test()`.
#' @references
#' DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or
#' more correlated receiver operating characteristic curves: a
#' nonparametric approach. *Biometrics*. 1988;44(3):837-845.
#' @examples
#' set.seed(1)
#' y <- rbinom(200, 1, 0.4)
#' d <- data.frame(y = y, m1 = y + rnorm(200), m2 = 0.5 * y + rnorm(200))
#' testroccompare(d, m1, m2, y, plot = FALSE)
#' @export
testroccompare <- function(data, predictor1, predictor2, outcome, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  p1 <- rlang::as_name(rlang::ensym(predictor1))
  p2 <- rlang::as_name(rlang::ensym(predictor2))
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  df <- drop_missing(data, c(p1, p2, outcome_nm), na.rm = na.rm)
  if (!is.numeric(df[[p1]]) || !is.numeric(df[[p2]])) stop("Both predictors must be numeric.", call. = FALSE)
  outcome_levels <- assert_two_groups(df, outcome_nm)
  df[[outcome_nm]] <- factor(df[[outcome_nm]], levels = outcome_levels)
  is_pos <- df[[outcome_nm]] == outcome_levels[2]
  n1 <- sum(is_pos); n0 <- sum(!is_pos)
  if (n1 < 2 || n0 < 2) stop("DeLong comparison needs at least two observations per outcome class.", call. = FALSE)

  placements <- function(x) {
    pos <- x[is_pos]; neg <- x[!is_pos]
    r_all <- rank(c(pos, neg))
    v10 <- (r_all[seq_len(n1)] - rank(pos)) / n0
    v01 <- 1 - (r_all[n1 + seq_len(n0)] - rank(neg)) / n1
    list(auc = mean(v10), v10 = v10, v01 = v01)
  }
  a <- placements(df[[p1]])
  b <- placements(df[[p2]])
  s10 <- stats::cov(cbind(a$v10, b$v10))
  s01 <- stats::cov(cbind(a$v01, b$v01))
  vmat <- s10 / n1 + s01 / n0
  diff <- a$auc - b$auc
  se_diff <- sqrt(max(0, vmat[1, 1] + vmat[2, 2] - 2 * vmat[1, 2]))
  if (se_diff == 0) stop("The two markers give identical rankings; the AUC difference has zero variance.", call. = FALSE)
  z <- stats::qnorm(1 - alpha / 2)
  z_stat <- diff / se_diff
  p_value <- 2 * stats::pnorm(-abs(z_stat))

  se_auc <- sqrt(diag(vmat))
  auc_table <- .bt(
    predictor = c(p1, p2),
    auc = c(a$auc, b$auc),
    se = se_auc,
    conf.low = c(a$auc, b$auc) - z * se_auc,
    conf.high = c(a$auc, b$auc) + z * se_auc,
    correlation = vmat[1, 2] / prod(se_auc)
  )
  primary <- .bt(
    method = "Paired DeLong test (AUC1 - AUC2)",
    estimate = diff,
    statistic = z_stat,
    p.value = p_value,
    conf.low = diff - z * se_diff,
    conf.high = diff + z * se_diff
  )
  h0 <- paste0("H0: AUC(", p1, ") = AUC(", p2, ") for ", outcome_nm, ".")
  effect <- .bt(name = paste0("AUC difference (", p1, " - ", p2, ")"), estimate = diff, magnitude = NA_character_)

  plt <- if (plot) {
    curve <- function(x, nm) {
      thr <- sort(unique(x), decreasing = TRUE)
      data.frame(
        predictor = nm,
        fpr = c(0, vapply(thr, function(t) mean(x[!is_pos] >= t), numeric(1))),
        tpr = c(0, vapply(thr, function(t) mean(x[is_pos] >= t), numeric(1)))
      )
    }
    cd <- rbind(curve(df[[p1]], p1), curve(df[[p2]], p2))
    lb <- stats::setNames(paste0(c(p1, p2), ", AUC ", format_stat(auc_table$auc), " (", format_stat(auc_table$conf.low), " to ",
                          format_stat(pmin(auc_table$conf.high, 1)), ")"), c(p1, p2))
    cd$predictor <- factor(cd$predictor, c(p1, p2))
    yd <- do.call(rbind, lapply(split(cd, cd$predictor), function(d) d[which.max(d$tpr - d$fpr), ]))
    ggplot2::ggplot(cd, ggplot2::aes(x = .data$fpr, y = .data$tpr, colour = .data$predictor)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
      ggplot2::geom_path(linewidth = 0.8) +
      ggplot2::geom_point(data = yd, size = 2.5, shape = 21, fill = "white", stroke = 1.2) +
      ggplot2::scale_colour_manual(values = stats::setNames(unname(TFCOL[c("a1", "a2")]), c(p1, p2)), labels = lb) +
      ggplot2::coord_equal() +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1)) +
      ggplot2::labs(subtitle = paste0("Paired DeLong test, difference ", format_stat(diff), ", ", pstr(p_value)),
                    caption = paste0("Open circles: Youden-optimal cut-offs. ", 100 * (1 - alpha), "% CIs."),
                    x = "1 - Specificity", y = "Sensitivity", colour = NULL) +
      tftheme() + ggplot2::theme(legend.position = "bottom")
  } else NULL

  out <- new_testflow(
    "roccompare",
    "paired ROC comparison",
    paste(p1, "vs", p2),
    outcome_nm,
    data = df,
    descriptives = descriptives_numeric(df, c(p1, p2), outcome_nm),
    assumptions = assumption_checks(check_independence_note("Both markers are measured on the same independent subjects; higher values indicate the positive class.")),
    recommended = list(test = "Paired DeLong test", rationale = "Two markers on the same subjects give correlated AUCs; the DeLong covariance accounts for that correlation."),
    primary_test = add_null_hypothesis(primary, h0),
    alternative_tests = list(auc_table = auc_table),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "roccompare"
  )
  out$interpretation <- make_report(out, alpha)
  out
}
