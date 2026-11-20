#' Test association between two categorical variables
#' @param formula A formula such as `x ~ y`, or a data frame when using pipe/data-first style.
#' @param data A data frame, or a first categorical column when using data-first style.
#' @param y Second categorical column. Optional when using formula style.
#' @param alpha Significance level.
#' @param fisher_threshold Expected-count threshold for Fisher's exact test.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param biascorrect Logical. If `FALSE` (default), the effect size is the
#'   usual Cramer's V, \eqn{V = \sqrt{\chi^2 / (n(\min(r, c) - 1))}}, which is
#'   biased upwards in small samples and sparse tables. If `TRUE`, the
#'   Bergsma (2013) bias-corrected V is reported instead.
#' @param event For 2 x 2 tables, the outcome level counted as the event.
#'   Default `NULL` uses the second level of the outcome.
#' @return A `testflow` object with class `testflow_categorical`. The object is
#' a list containing the cleaned data, categorical descriptives, assumption
#' checks, recommended association test, primary test result with null
#' hypothesis, alternative chi-square and Fisher results, effect size, optional
#' `ggplot`, original call, and report text.
#'
#' `$posthoc` holds one row per cell: `observed`, `expected`, the adjusted
#' standardised residual `stdres` (as `stats::chisq.test()$stdres`), its
#' two-sided normal `p.value` and the Benjamini-Hochberg `p.adj`.
#'
#' For a 2 x 2 table, `$alternative_tests$effect_measures` compares the
#' second level of the group (columns) with the first, for the `event` level
#' of the outcome (rows): risk difference `RD` with the Newcombe hybrid score
#' interval (Wilson limits of each risk), risk ratio `RR` (log scale,
#' Katz), odds ratio `OR` (log scale, Woolf) and the number needed to treat
#' `NNT` = 1/|RD| with limits 1/|RD limits| (Altman 1998). When the RD
#' interval includes 0, the NNT interval runs through infinity and `note`
#' says so. With a zero cell, 0.5 is added to every cell for RR and OR.
#' @references
#' Pearson, K. (1900). On the criterion that a given system of deviations from
#' the probable in the case of a correlated system of variables is such that it
#' can be reasonably supposed to have arisen from random sampling.
#' \emph{Philosophical Magazine}, 50(302), 157-175.
#'
#' Fisher, R. A. (1922). On the interpretation of chi-square from contingency
#' tables, and the calculation of P. \emph{Journal of the Royal Statistical
#' Society}, 85(1), 87-94.
#'
#' Cramer, H. (1946). \emph{Mathematical Methods of Statistics}. Princeton.
#'
#' Newcombe, R. G. (1998). Interval estimation for the difference between
#' independent proportions: comparison of eleven methods. \emph{Statistics
#' in Medicine}, 17(8), 873-890.
#'
#' Altman, D. G. (1998). Confidence intervals for the number needed to treat.
#' \emph{BMJ}, 317(7168), 1309-1312.
#'
#' Bergsma, W. (2013). A bias-correction for Cramer's V and Tschuprow's T.
#' \emph{Journal of the Korean Statistical Society}, 42(3), 323-328.
#' @export
testcategorical <- function(formula, data, y = NULL, alpha = 0.05, fisher_threshold = 5, plot = TRUE, na.rm = TRUE, biascorrect = FALSE,
                            event = NULL) {
  vars <- resolve_formula_pair(tf_formula_arg(substitute(formula), parent.frame()), tf_formula_arg(substitute(data), parent.frame()), substitute(y), missing(y), rhs_label = "y")
  x_nm <- vars$outcome
  y_nm <- vars$group
  data_obj <- resolve_data_first_or_formula(formula, data)
  df <- drop_missing(data_obj, c(x_nm, y_nm), na.rm = na.rm)
  tab <- table(df[[x_nm]], df[[y_nm]])
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  fisher <- .fisher_safe(tab)
  use_fisher <- any(chi$expected < fisher_threshold)
  fisher_label <- if (grepl("Monte Carlo", fisher$method, fixed = TRUE)) fisher$method else "Fisher exact test"
  recommendation <- if (use_fisher) fisher_label else "Chi-square test of independence"
  primary <- if (use_fisher) fisher else chi
  effect <- cramers_v(tab, biascorrect = biascorrect)
  expected <- check_expected_counts(tab, fisher_threshold)
  plt <- if (plot) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_nm]], fill = .data[[y_nm]])) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(title = "Categorical association workflow", subtitle = plot_subtitle(recommendation, primary), caption = paste0(effect$name[1], " = ", format_stat(effect$estimate[1]), ", ", effect$magnitude[1]), x = x_nm, y = "Proportion", fill = y_nm) +
      ggplot2::theme_minimal()
  } else NULL
  h0 <- h0_no_association(x_nm, y_nm)
  out <- new_testflow("categorical", "two categorical variables", x_nm, y_nm, data = df, descriptives = descriptives_categorical(df, c(x_nm, y_nm)), assumptions = assumption_checks(check_independence_note(), expected), recommended = list(test = recommendation), primary_test = add_null_hypothesis(safe_tidy_htest(primary, recommendation), h0), alternative_tests = list(chi_square = add_null_hypothesis(safe_tidy_htest(chi, "Chi-square test of independence"), h0), fisher = add_null_hypothesis(safe_tidy_htest(fisher, fisher_label), h0),
                      effect_measures = if (all(dim(tab) == 2)) effect_measures_2x2(tab, event, alpha)), posthoc = cell_residuals(chi),
                      effect_size = effect, plot = plt, call = match.call(), subclass = "categorical")
  out$interpretation <- make_report(out, alpha)
  out
}

# Cell-wise adjusted standardised residuals (Agresti 2013, Sec 2.4.5).
cell_residuals <- function(chi) {
  o <- as.data.frame(chi$observed, stringsAsFactors = FALSE); z <- as.vector(chi$stdres)
  p <- 2 * stats::pnorm(-abs(z))
  .bt(row = o[[1]], column = o[[2]], observed = o$Freq, expected = as.vector(chi$expected), stdres = z,
      p.value = p, p.adj = stats::p.adjust(p, "BH"), p.adjust.method = "BH")
}

# RD (Newcombe), RR (Katz), OR (Woolf) and NNT (Altman) for group column 2 vs 1.
effect_measures_2x2 <- function(tab, event = NULL, alpha = 0.05) {
  ev <- event %||% rownames(tab)[2]
  if (!ev %in% rownames(tab)) stop("`event` must be a level of the outcome.", call. = FALSE)
  e <- tab[ev, ]; n <- colSums(tab); p <- e / n; z <- stats::qnorm(1 - alpha / 2)
  w1 <- wilson_ci(e[[1]], n[[1]], alpha); w2 <- wilson_ci(e[[2]], n[[2]], alpha)
  rd <- p[[2]] - p[[1]]
  rdci <- c(rd - sqrt((p[[2]] - w2[1])^2 + (w1[2] - p[[1]])^2), rd + sqrt((w2[2] - p[[2]])^2 + (p[[1]] - w1[1])^2))
  cc <- if (any(tab == 0)) 0.5 else 0
  a <- e + cc; m <- n + 2 * cc; b <- m - a
  rr <- (a[[2]] / m[[2]]) / (a[[1]] / m[[1]]); se_rr <- sqrt(sum(1 / a - 1 / m))
  or <- (a[[2]] / b[[2]]) / (a[[1]] / b[[1]]); se_or <- sqrt(sum(1 / a + 1 / b))
  spans <- rdci[1] <= 0 && rdci[2] >= 0
  nnt_ci <- if (spans) c(1 / abs(rdci[1]), 1 / abs(rdci[2])) else sort(1 / abs(rdci))
  dir <- if (rd < 0) "benefit (fewer events)" else "harm (more events)"
  cmp <- paste0(colnames(tab)[2], " vs ", colnames(tab)[1])
  .bt(measure = c("RD", "RR", "OR", "NNT"), estimate = c(rd, rr, or, 1 / abs(rd)),
      conf.low = c(rdci[1], exp(log(rr) - z * se_rr), exp(log(or) - z * se_or), nnt_ci[1]),
      conf.high = c(rdci[2], exp(log(rr) + z * se_rr), exp(log(or) + z * se_or), nnt_ci[2]),
      method = c("Newcombe hybrid score", "Log (Katz)", "Log (Woolf)", "1/|RD| (Altman)"),
      comparison = cmp, event = ev,
      note = c(NA, if (cc > 0) rep("0.5 added to each cell (zero cell)", 2) else rep(NA, 2),
               if (spans) paste0("RD CI includes 0: NNT runs from ", format_stat(nnt_ci[1]), " (", if (rd < 0) "NNTB" else "NNTH",
                                 ") through infinity to ", format_stat(nnt_ci[2]), " (", if (rd < 0) "NNTH" else "NNTB", ")")
               else paste0("NNT for ", dir)))
}
