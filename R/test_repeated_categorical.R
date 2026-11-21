#' Test repeated categorical measurements
#' @param data A data frame.
#' @param measures Repeated binary columns selected with tidyselect syntax.
#' @param id Optional subject identifier.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param event Level whose proportion is plotted. Default `NULL` uses the
#'   second level over all measures (the level coded 1 for the test).
#'   Cochran's Q does not depend on this choice; only the plot changes.
#' @return A `testflow` object with class `testflow_repeated_categorical`. The
#' object is a list containing the cleaned repeated categorical data,
#' descriptive counts, assumption checks, Cochran Q test result with null
#' hypothesis, pairwise McNemar post-hoc comparisons, effect size, optional
#' `ggplot`, original call, and report text.
#' @references
#' Cochran, W. G. (1950). The comparison of percentages in matched samples.
#' \emph{Biometrika}, 37(3/4), 256-266.
#'
#' McNemar, Q. (1947). Note on the sampling error of the difference between
#' correlated proportions or percentages. \emph{Psychometrika}, 12(2), 153-157.
#' @export
testrepeatedcategorical <- function(data, measures, id = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE, event = NULL) {
  measure_nms <- tidyselect_names(data, {{ measures }})
  warn_if_two_measures_repeated(measure_nms)
  for (nm in measure_nms) warn_if_nonbinary(data[[nm]], paste0("`", nm, "`"))
  df <- drop_missing(data, measure_nms, na.rm = na.rm)
  mat <- as.matrix(df[, measure_nms, drop = FALSE])
  mat <- apply(mat, 2, function(x) as.integer(as.factor(x)) - 1L)
  warn_if(!all(mat %in% c(0, 1)), "Repeated categorical measures should be binary after recoding.")
  test <- cochran_q_test(mat)
  col_sums <- colSums(mat)
  counts <- .bt(time = names(col_sums), success = as.integer(col_sums))
  effect <- cochran_q_effect_size(test, n = nrow(mat), k = ncol(mat))
  plt <- if (plot) {
    v <- unlist(lapply(df[measure_nms], function(x) as.character(catlev(x))))
    lv <- levels(catlev(unlist(lapply(df[measure_nms], catlev)))); ev <- event %||% lv[min(2, length(lv))]
    if (is.logical(ev)) ev <- c("No", "Yes")[ev + 1]
    k <- vapply(df[measure_nms], function(x) sum(as.character(catlev(x)) == ev), 1); n <- nrow(df)
    ci <- t(vapply(k, function(x) wilson_ci(x, n, alpha), numeric(2)))
    pd <- data.frame(time = factor(measure_nms, measure_nms), p = k / n, lo = ci[, 1], hi = ci[, 2])
    ggplot2::ggplot(pd, ggplot2::aes(x = .data$time, y = .data$p, group = 1)) +
      ggplot2::geom_line(colour = TFCOL[["a1"]]) +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), colour = TFCOL[["a1"]]) +
      ggplot2::scale_y_continuous(labels = function(b) paste0(100 * b, "%")) +
      ggplot2::labs(subtitle = plot_subtitle("Cochran Q test", test), caption = paste0(
        "Points: percent with ", ev, " and ", 100 * (1 - alpha), "% Wilson CI (n = ", n, ")."),
        x = NULL, y = paste0("Percent ", ev)) +
      tftheme()
  } else NULL
  h0 <- paste0("H0: the success proportions are equal across repeated categorical measures.")
  out <- new_testflow("repeated_categorical", "repeated categorical measurements", paste(measure_nms, collapse = ", "), data = df, descriptives = counts, assumptions = assumption_checks(assumption_check("Repeated binary measurements", "assumed", "Same subjects should be measured at 3 or more time points."), assumption_check("Complete repeated data", ifelse(nrow(df) == nrow(mat), "acceptable", "warning"), "Missingness should be handled explicitly or via complete-case analysis.")), recommended = list(test = "Cochran Q test"), primary_test = add_null_hypothesis(safe_tidy_htest(test, "Cochran Q test"), h0), alternative_tests = list(pairwise_mcnemar = pairwise_mcnemar(mat, measure_nms)), posthoc = pairwise_mcnemar(mat, measure_nms), effect_size = effect, plot = plt, call = match.call(), subclass = "repeated_categorical")
  out$interpretation <- make_report(out, alpha)
  out
}
