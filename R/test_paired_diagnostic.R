#' Compare the accuracy of two binary tests on the same subjects
#'
#' Sensitivity and specificity of two tests applied to the same subjects,
#' each compared with McNemar's test within the diseased (sensitivity) and
#' non-diseased (specificity) subjects.
#'
#' @param data A data frame.
#' @param test1,test2 Binary test result columns; the second level of each is
#'   the positive result.
#' @param reference Binary reference standard column; the second level is
#'   disease present.
#' @param exact Logical; use the exact binomial McNemar test (conditional on
#'   the discordant pairs) instead of the chi-square version.
#' @param correct Logical; continuity correction of the chi-square McNemar
#'   test (as `stats::mcnemar.test()`).
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_paired_diagnostic`.
#'   `$primary_test` has two rows (`measure` = `"Sensitivity"`,
#'   `"Specificity"`) with the difference test 2 minus test 1 (`estimate`),
#'   its Wald interval for paired proportions (`conf.low`, `conf.high`), and
#'   the McNemar `statistic`, `parameter` and `p.value`.
#'   `$alternative_tests$comparison` adds each test's value (`test1`,
#'   `test2`), the subjects in that disease group (`n`) and the discordant
#'   counts (`only_test1`: test 1 correct and test 2 wrong; `only_test2`:
#'   the reverse). `$alternative_tests$accuracy` gives Se and Sp of each test
#'   with Clopper-Pearson intervals.
#' @details
#' Among the \eqn{n} subjects of one disease group, let \eqn{b} and \eqn{c}
#' be the pairs where only test 1 or only test 2 is correct. The difference
#' is \eqn{(c - b)/n} with standard error
#' \eqn{\sqrt{b + c - (c - b)^2/n}/n}, and McNemar's statistic is
#' \eqn{(|b - c| - 1)^2/(b + c)} (without the 1 when `correct = FALSE`), or
#' the exact binomial test of \eqn{b} out of \eqn{b + c} with \eqn{p = 0.5}.
#' @references
#' Hawass NED. Comparing the sensitivities and specificities of two
#' diagnostic procedures performed on the same group of patients. *British
#' Journal of Radiology*. 1997;70:360-366.
#'
#' Pepe MS. *The Statistical Evaluation of Medical Tests for Classification
#' and Prediction*. Oxford University Press; 2003.
#' @examples
#' set.seed(1)
#' d <- data.frame(dis = rbinom(200, 1, 0.4))
#' d$t1 <- ifelse(d$dis == 1, rbinom(200, 1, 0.8), rbinom(200, 1, 0.1))
#' d$t2 <- ifelse(d$dis == 1, rbinom(200, 1, 0.9), rbinom(200, 1, 0.2))
#' testpaireddiagnostic(d, t1, t2, dis, plot = FALSE)
#' @export
testpaireddiagnostic <- function(data, test1, test2, reference, exact = FALSE, correct = TRUE, alpha = 0.05,
                                 plot = TRUE, na.rm = TRUE) {
  t1 <- rlang::as_name(rlang::ensym(test1)); t2 <- rlang::as_name(rlang::ensym(test2))
  rf <- rlang::as_name(rlang::ensym(reference))
  df <- drop_missing(data, c(t1, t2, rf), na.rm = na.rm)
  pos <- function(nm) df[[nm]] == assert_two_groups(df, nm)[2]
  a <- pos(t1); b <- pos(t2); d <- pos(rf); z <- stats::qnorm(1 - alpha / 2)
  one <- function(grp, ok1, ok2, measure) {
    n <- sum(grp); x1 <- sum(ok1 & !ok2 & grp); x2 <- sum(!ok1 & ok2 & grp)
    est <- (x2 - x1) / n; se <- sqrt(x1 + x2 - (x2 - x1)^2 / n) / n
    if (exact) {
      tt <- stats::binom.test(x1, x1 + x2); st <- NA_real_; par <- NA_real_
      m <- "Exact McNemar test (binomial)"
    } else {
      st <- (abs(x1 - x2) - if (correct) 1 else 0)^2 / (x1 + x2); par <- 1
      tt <- list(p.value = stats::pchisq(st, 1, lower.tail = FALSE))
      m <- if (correct) "McNemar test with continuity correction" else "McNemar test"
    }
    .bt(measure = measure, method = m, test1 = sum(ok1 & grp) / n, test2 = sum(ok2 & grp) / n, estimate = est,
        conf.low = est - z * se, conf.high = est + z * se, statistic = st, parameter = par, p.value = tt$p.value,
        n = n, only_test1 = x1, only_test2 = x2)
  }
  comparison <- rbind(one(d, a, b, "Sensitivity"), one(!d, !a, !b, "Specificity"))
  cp <- function(x, m) as.numeric(stats::binom.test(x, m, conf.level = 1 - alpha)$conf.int)
  acc <- .map_dfr(list(c(t1, "a"), c(t2, "b")), function(p) {
    r <- if (p[2] == "a") a else b
    se <- cp(sum(r & d), sum(d)); sp <- cp(sum(!r & !d), sum(!d))
    .bt(test = p[1], measure = c("Sensitivity", "Specificity"), estimate = c(mean(r[d]), mean(!r[!d])),
        conf.low = c(se[1], sp[1]), conf.high = c(se[2], sp[2]))
  })
  primary <- comparison[, c("measure", "method", "estimate", "conf.low", "conf.high", "statistic", "parameter", "p.value")]
  primary$conf.level <- 1 - alpha
  h0 <- paste0("H0: ", t1, " and ", t2, " have equal sensitivity (diseased) and equal specificity (non-diseased).")
  plt <- if (plot) {
    ggplot2::ggplot(acc, ggplot2::aes(x = .data$test, y = .data$estimate, colour = .data$test)) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high)) +
      ggplot2::facet_wrap(~measure) + ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = "Paired diagnostic accuracy workflow", subtitle = paste0("McNemar p: Se ",
                    format_p(comparison$p.value[1]), ", Sp ", format_p(comparison$p.value[2])), x = NULL, y = NULL) +
      ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none")
  } else NULL
  out <- new_testflow("paired_diagnostic", "two binary tests on the same subjects", paste(t1, "vs", t2), rf, data = df,
                      descriptives = descriptives_categorical(df, c(t1, t2, rf)),
                      assumptions = assumption_checks(check_independence_note(
                        "Subjects are independent; both tests and the reference are read on every subject.")),
                      recommended = list(test = comparison$method[1], rationale =
                        "Paired binary results: discordant pairs within diseased and non-diseased subjects drive the comparison."),
                      primary_test = add_null_hypothesis(primary, h0),
                      alternative_tests = list(comparison = comparison, accuracy = acc),
                      effect_size = .bt(name = c("Sensitivity difference", "Specificity difference"), estimate = comparison$estimate,
                                        magnitude = NA_character_),
                      plot = plt, call = match.call(), subclass = "paired_diagnostic")
  pct <- function(x) paste0(format_stat(100 * x, 1), "%")
  txt <- function(i) paste0(comparison$measure[i], " was ", pct(comparison$test1[i]), " for ", t1, " and ", pct(comparison$test2[i]),
                            " for ", t2, " (difference ", pct(comparison$estimate[i]), ", ", primary_ci_level_text(out), " CI ",
                            pct(comparison$conf.low[i]), " to ", pct(comparison$conf.high[i]), "; ", comparison$method[i], " p = ",
                            format_p(comparison$p.value[i]), ").")
  out$interpretation <- paste(txt(1), txt(2))
  out
}
