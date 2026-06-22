#' Return a ready-to-use testflow report
#' @param x A testflow object.
#' @param ... Unused.
#' @export
report <- function(x, ...) {
  UseMethod("report")
}

#' @export
report.testflow <- function(x, ...) {
  x$interpretation %||% make_report(x)
}

#' Return a ready-to-use testflow report
#' @param x A testflow object.
#' @export
report_test <- function(x) {
  report(x)
}

make_report <- function(x, alpha = 0.05) {
  if (identical(x$workflow, "correlation_matrix")) {
    effect <- x$effect_size
    effect_text <- if (is.null(effect) || nrow(effect) == 0 || is.na(effect$estimate[1]) || is.na(effect$magnitude[1])) {
      ""
    } else {
      paste0(" The strongest absolute correlation was ", format_stat(effect$estimate[1]), " (", effect$magnitude[1], ").")
    }
    return(paste0(
      "The correlation matrix workflow for ", x$outcome %||% "the selected variables",
      " reported pairwise correlations using ", x$recommended$test %||% x$recommended,
      "; the smallest pairwise p-value was ", format_p(primary_p(x)), ".", effect_text
    ))
  }

  if (identical(x$workflow, "outliers")) {
    n <- primary_statistic(x)
    return(paste0("The outlier workflow flagged ", format_count(n), " ", ifelse(n == 1, "row", "rows"), " for review."))
  }

  p <- primary_p(x)
  method <- x$recommended$test %||% x$recommended
  stat <- primary_statistic(x)
  df <- primary_df(x)
  effect <- x$effect_size
  effect_text <- if (is.null(effect) || nrow(effect) == 0 || is.na(effect$estimate[1]) || is.na(effect$magnitude[1])) {
    ""
  } else {
    paste0(" The effect size was ", effect$magnitude[1], " (", effect$name[1], " = ", format_stat(effect$estimate[1]), ").")
  }
  ci <- primary_ci(x)
  ci_text <- if (any(is.na(ci))) "" else paste0(" The 95% confidence interval was [", format_stat(ci[1]), ", ", format_stat(ci[2]), "].")
  h0 <- primary_h0(x)
  h0_text <- if (is.na(h0)) "" else paste0(" ", h0)
  result <- ifelse(p < alpha, "showed a statistically significant result", "did not show a statistically significant result")
  paste0(
    "The ", x$design, " workflow for ", x$outcome %||% "the outcome",
    " ", result, " using ", method, ", statistic = ", format_stat(stat),
    ifelse(is.na(df), "", paste0(", df = ", format_stat(df))),
    ", p = ", format_p(p), ".", ci_text, effect_text, h0_text
  )
}
