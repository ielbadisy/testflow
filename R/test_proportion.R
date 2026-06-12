#' Test a one-sample proportion
#' @param data A data frame.
#' @param outcome Categorical outcome column.
#' @param success Value counted as success.
#' @param p Reference probability.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_proportion <- function(data, outcome, success, p = 0.5, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  df <- drop_missing(data, outcome_nm, na.rm = na.rm)
  success_n <- sum(df[[outcome_nm]] == success, na.rm = TRUE)
  total_n <- sum(!is.na(df[[outcome_nm]]))
  binom <- stats::binom.test(success_n, total_n, p = p)
  prop <- stats::prop.test(success_n, total_n, p = p, correct = FALSE)
  effect <- tibble::tibble(name = "Observed proportion", estimate = success_n / total_n, magnitude = NA_character_)
  plt <- if (plot) {
    counts <- c(success_n, total_n - success_n)
    pd <- tibble::tibble(level = c("success", "other"), n = counts, prop = counts / total_n)
    ggplot2::ggplot(pd, ggplot2::aes(x = .data$level, y = .data$prop, fill = .data$level)) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept = p, linetype = "dashed") +
      ggplot2::labs(title = "One-sample proportion workflow", subtitle = plot_subtitle("Exact binomial test", binom), x = NULL, y = "Proportion") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  } else NULL
  h0 <- h0_proportion(outcome_nm, p)
  out <- new_testflow("proportion", "one categorical proportion", outcome_nm, data = df, descriptives = descriptives_categorical(df, outcome_nm), recommended = list(test = "Exact binomial test"), primary_test = add_null_hypothesis(safe_tidy_htest(binom, "Exact binomial test"), h0), alternative_tests = list(prop_test = add_null_hypothesis(safe_tidy_htest(prop, "One-sample proportion test"), h0)), effect_size = effect, plot = plt, call = match.call(), subclass = "proportion")
  out$interpretation <- make_report(out, alpha)
  out
}
