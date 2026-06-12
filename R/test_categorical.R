#' Test association between two categorical variables
#' @param data A data frame.
#' @param x First categorical column.
#' @param y Second categorical column.
#' @param alpha Significance level.
#' @param fisher_threshold Expected-count threshold for Fisher's exact test.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_categorical <- function(data, x, y, alpha = 0.05, fisher_threshold = 5, plot = TRUE, na.rm = TRUE) {
  x_nm <- rlang::as_name(rlang::ensym(x)); y_nm <- rlang::as_name(rlang::ensym(y))
  df <- drop_missing(data, c(x_nm, y_nm), na.rm = na.rm)
  tab <- table(df[[x_nm]], df[[y_nm]])
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  fisher <- stats::fisher.test(tab)
  recommendation <- if (any(chi$expected < fisher_threshold)) "Fisher exact test" else "Chi-square test of independence"
  primary <- if (recommendation == "Fisher exact test") fisher else chi
  effect <- cramers_v(tab)
  plt <- if (plot) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_nm]], fill = .data[[y_nm]])) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(title = "Categorical association workflow", subtitle = plot_subtitle(recommendation, primary), caption = paste0(effect$name[1], " = ", format_stat(effect$estimate[1]), ", ", effect$magnitude[1]), x = x_nm, y = "Proportion", fill = y_nm) +
      ggplot2::theme_minimal()
  } else NULL
  out <- new_testflow("categorical", "two categorical variables", x_nm, y_nm, data = df, descriptives = descriptives_categorical(df, c(x_nm, y_nm)), assumptions = list("Expected counts" = tibble::as_tibble(chi$expected, rownames = x_nm)), recommended = list(test = recommendation), primary_test = safe_tidy_htest(primary, recommendation), alternative_tests = list(chi_square = safe_tidy_htest(chi, "Chi-square test of independence"), fisher = safe_tidy_htest(fisher, "Fisher exact test")), effect_size = effect, plot = plt, call = match.call(), subclass = "categorical")
  out$interpretation <- make_report(out, alpha)
  out
}
