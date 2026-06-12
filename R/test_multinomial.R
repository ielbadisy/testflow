#' Test multinomial goodness of fit
#' @param data A data frame.
#' @param outcome Categorical outcome column.
#' @param p Expected probabilities, or `NULL` for equal probabilities.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_multinomial <- function(data, outcome, p = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  df <- drop_missing(data, outcome_nm, na.rm = na.rm)
  counts <- as.numeric(table(df[[outcome_nm]]))
  levels <- names(table(df[[outcome_nm]]))
  if (is.null(p)) p <- rep(1 / length(counts), length(counts))
  expected <- sum(counts) * p
  chisq <- suppressWarnings(stats::chisq.test(x = counts, p = p))
  pairwise <- purrr::map_dfr(seq_along(counts), function(i) {
    bt <- stats::binom.test(counts[i], sum(counts), p = p[i])
    tibble::tibble(level = levels[i], observed = counts[i], expected = expected[i], p = bt$p.value)
  })
  effect <- tibble::tibble(name = "Chi-square goodness-of-fit", estimate = unname(chisq$statistic), magnitude = NA_character_)
  plt <- if (plot) {
    pd <- tibble::tibble(level = levels, observed = counts, expected = expected) |>
      tidyr::pivot_longer(c("observed", "expected"), names_to = "type", values_to = "n")
    ggplot2::ggplot(pd, ggplot2::aes(x = .data$level, y = .data$n, fill = .data$type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::labs(title = "Multinomial workflow", subtitle = plot_subtitle("Chi-square goodness-of-fit", chisq), x = outcome_nm, y = "Count") +
      ggplot2::theme_minimal()
  } else NULL
  out <- new_testflow("multinomial", "one multinomial categorical variable", outcome_nm, data = df, descriptives = descriptives_categorical(df, outcome_nm), recommended = list(test = "Chi-square goodness-of-fit"), primary_test = safe_tidy_htest(chisq, "Chi-square goodness-of-fit"), alternative_tests = list(pairwise_binomial = pairwise), effect_size = effect, plot = plt, call = match.call(), subclass = "multinomial")
  out$posthoc <- pairwise
  out$interpretation <- make_report(out, alpha)
  out
}
