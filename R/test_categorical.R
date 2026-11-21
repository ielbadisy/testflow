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
  em <- if (all(dim(tab) == 2)) effect_measures_2x2(tab, event, alpha)
  plt <- if (plot) stack_plot(df[[y_nm]], df[[x_nm]], y_nm, x_nm, plot_subtitle(recommendation, primary),
                              paste(c(effstr(effect), emstr(em, alpha)), collapse = "\n"), event = em$event[1]) else NULL
  h0 <- h0_no_association(x_nm, y_nm)
  out <- new_testflow("categorical", "two categorical variables", x_nm, y_nm, data = df, descriptives = descriptives_categorical(df, c(x_nm, y_nm)), assumptions = assumption_checks(check_independence_note(), expected), recommended = list(test = recommendation), primary_test = add_null_hypothesis(safe_tidy_htest(primary, recommendation), h0), alternative_tests = list(chi_square = add_null_hypothesis(safe_tidy_htest(chi, "Chi-square test of independence"), h0), fisher = add_null_hypothesis(safe_tidy_htest(fisher, fisher_label), h0),
                      effect_measures = em), posthoc = cell_residuals(chi),
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

# Logical and 0/1 levels read as words; other levels keep their order.
catlev <- function(x) {
  if (is.logical(x)) return(factor(x, c(FALSE, TRUE), c("No", "Yes")))
  factor(x, lvls(x))
}

# Bars per group (x) stacked by outcome (fill), proportions within each group, with the count in each segment.
stack_plot <- function(g, y, gname, yname, subtitle = NULL, caption = NULL, event = NULL) {
  if (is.logical(y) && !is.null(event)) event <- c("FALSE" = "No", "TRUE" = "Yes")[[as.character(event)]]
  g <- catlev(g); y <- catlev(y)
  # The event level goes at the bottom of each bar, in blue, so its proportion reads off the axis.
  lv <- levels(y); if (!is.null(event) && event %in% lv) lv <- c(setdiff(lv, event), event)
  y <- factor(y, lv); k <- length(lv)
  tb <- as.data.frame(table(g = g, y = y)); tb$p <- stats::ave(tb$Freq, tb$g, FUN = function(f) f / max(sum(f), 1))
  tb$ink <- ifelse(as.integer(tb$y) > k / 2, "white", "grey15")
  ggplot2::ggplot(tb, ggplot2::aes(x = .data$g, y = .data$p, fill = .data$y)) +
    ggplot2::geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data$p >= 0.06, .data$Freq, ""), colour = .data$ink),
                       position = ggplot2::position_stack(vjust = 0.5), size = 3) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_manual(values = stats::setNames(tfgrey(k), lv)) +
    ggplot2::scale_x_discrete(labels = wraplab(14)) +
    ggplot2::scale_y_continuous(labels = function(b) paste0(100 * b, "%"), expand = ggplot2::expansion(c(0, 0.02))) +
    ggplot2::labs(subtitle = subtitle, caption = caption, x = gname, y = paste0("Percent within ", gname), fill = yname) +
    tftheme() + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}

# One line with RD, RR and OR (and their CIs) of a 2 x 2 effect-measure table.
emstr <- function(em, alpha = 0.05) {
  if (is.null(em)) return(NULL)
  f <- function(i, pct = FALSE) { k <- if (pct) 100 else 1
    paste0(em$measure[i], " ", format_stat(k * em$estimate[i], if (pct) 1 else 2), " (", format_stat(k * em$conf.low[i], if (pct) 1 else 2),
           " to ", format_stat(k * em$conf.high[i], if (pct) 1 else 2), ")") }
  paste0(em$comparison[1], ", event ", em$event[1], " (", 100 * (1 - alpha), "% CI):\n", f(1, TRUE), " points; ", f(2),
         "; ", f(3))
}

#' Plot a categorical association test
#'
#' @param x A `testflow_categorical` object from [testcategorical()].
#' @param type `"bars"` (default) draws the stored plot (percent within each
#'   group); `"residuals"` draws the adjusted standardised residual of every
#'   cell (`$posthoc$stdres`), blue below and vermillion above expectation,
#'   with cells beyond +-1.96 in bold. Useful for r x c tables.
#' @param ... Passed to [plot.testflow()] (`title`, `subtitle`, `caption`).
#' @return A `ggplot` object.
#' @export
plot.testflow_categorical <- function(x, type = c("bars", "residuals"), ...) {
  if (match.arg(type) == "residuals") {
    r <- as.data.frame(x$posthoc); d <- x$data
    # Same level order as table() used for the test.
    r$row <- factor(r$row, rev(levels(as.factor(d[[x$outcome]])))); r$column <- factor(r$column, levels(as.factor(d[[x$group]])))
    lim <- max(2.5, abs(r$stdres)); sig <- abs(r$stdres) > stats::qnorm(0.975)
    x$plot <- ggplot2::ggplot(r, ggplot2::aes(x = .data$column, y = .data$row, fill = .data$stdres)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = formatC(.data$stdres, format = "f", digits = 1)), size = 3.3,
                         fontface = ifelse(sig, "bold", "plain"),
                         colour = ifelse(abs(r$stdres) > 0.6 * lim, "white", ifelse(sig, "black", "grey35"))) +
      ggplot2::scale_fill_gradient2(low = TFCOL[["a1"]], mid = "white", high = TFCOL[["a2"]], limits = c(-lim, lim)) +
      ggplot2::scale_x_discrete(labels = wraplab(14), position = "top") +
      ggplot2::labs(subtitle = plot_subtitle(x$recommended$test, x$primary_test), x = x$group, y = x$outcome,
                    fill = "Residual", caption = "Adjusted standardised residuals; bold: |z| > 1.96 (unadjusted).") +
      tftheme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  plot.testflow(x, ...)
}
