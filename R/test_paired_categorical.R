#' Test paired categorical measurements
#' @param data A data frame.
#' @param before Before categorical column.
#' @param after After categorical column.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_paired_categorical <- function(data, before, after, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  before_nm <- rlang::as_name(rlang::ensym(before)); after_nm <- rlang::as_name(rlang::ensym(after))
  df <- drop_missing(data, c(before_nm, after_nm), na.rm = na.rm)
  tab <- table(df[[before_nm]], df[[after_nm]])
  test <- stats::mcnemar.test(tab)
  disc <- tibble::tibble(before_only = if (all(dim(tab) == 2)) tab[2, 1] else NA_integer_, after_only = if (all(dim(tab) == 2)) tab[1, 2] else NA_integer_)
  effect <- tibble::tibble(name = "Discordant pairs", estimate = sum(disc, na.rm = TRUE), magnitude = NA_character_)
  plt <- if (plot) {
    long <- tidyr::pivot_longer(df, c(dplyr::all_of(before_nm), dplyr::all_of(after_nm)), names_to = "time", values_to = "value")
    ggplot2::ggplot(long, ggplot2::aes(x = .data$time, fill = .data$value)) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(title = "Paired categorical workflow", subtitle = plot_subtitle("McNemar test", test), x = NULL, y = "Proportion") +
      ggplot2::theme_minimal()
  } else NULL
  out <- new_testflow("paired_categorical", "paired categorical measurements", paste(before_nm, after_nm, sep = " -> "), data = df, descriptives = descriptives_categorical(df, c(before_nm, after_nm)), recommended = list(test = "McNemar test"), primary_test = safe_tidy_htest(test, "McNemar test"), alternative_tests = list(discordant_pairs = disc), effect_size = effect, plot = plt, call = match.call(), subclass = "paired_categorical")
  out$interpretation <- make_report(out, alpha)
  out
}
