#' Test paired categorical measurements
#' @param data A data frame.
#' @param before Before categorical column.
#' @param after After categorical column.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_paired_categorical`. The
#' object is a list containing the cleaned paired categorical data, categorical
#' descriptives, assumption checks, McNemar test result, discordant-pair table,
#' optional `ggplot`, original call, and report text.
#' @export
test_paired_categorical <- function(data, before, after, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  before_nm <- rlang::as_name(rlang::ensym(before)); after_nm <- rlang::as_name(rlang::ensym(after))
  df <- drop_missing(data, c(before_nm, after_nm), na.rm = na.rm)
  tab <- table(df[[before_nm]], df[[after_nm]])
  warn_if(!all(dim(tab) == 2), "McNemar's test requires a 2x2 table; the current input is not a binary paired design.")
  test <- stats::mcnemar.test(tab)
  discordant <- if (all(dim(tab) == 2)) sum(tab[c(1, 2), c(2, 1)]) else NA_integer_
  disc <- tibble::tibble(before_only = if (all(dim(tab) == 2)) tab[2, 1] else NA_integer_, after_only = if (all(dim(tab) == 2)) tab[1, 2] else NA_integer_)
  effect <- tibble::tibble(name = "Discordant pairs", estimate = sum(disc, na.rm = TRUE), magnitude = NA_character_)
  plt <- if (plot) {
    long <- tidyr::pivot_longer(df, c(dplyr::all_of(before_nm), dplyr::all_of(after_nm)), names_to = "time", values_to = "value")
    ggplot2::ggplot(long, ggplot2::aes(x = .data$time, fill = .data$value)) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(title = "Paired categorical workflow", subtitle = plot_subtitle("McNemar test", test), x = NULL, y = "Proportion") +
      ggplot2::theme_minimal()
  } else NULL
  discordant_message <- dplyr::case_when(
    is.na(discordant) ~ "Discordant pairs could not be counted.",
    discordant >= 10 ~ "Discordant pairs are adequate for the McNemar approximation.",
    TRUE ~ "Small numbers of discordant pairs favor exact McNemar."
  )
  out <- new_testflow("paired_categorical", "paired categorical measurements", paste(before_nm, after_nm, sep = " -> "), data = df, descriptives = descriptives_categorical(df, c(before_nm, after_nm)), assumptions = assumption_checks(assumption_check("Paired binary measurements", "assumed", "Same subjects should be measured twice."), assumption_check("Discordant pairs", ifelse(is.na(discordant) || discordant >= 10, "acceptable", "warning"), discordant_message, details = paste0("discordant = ", discordant))), recommended = list(test = "McNemar test"), primary_test = safe_tidy_htest(test, "McNemar test"), alternative_tests = list(discordant_pairs = disc), effect_size = effect, plot = plt, call = match.call(), subclass = "paired_categorical")
  out$interpretation <- make_report(out, alpha)
  out
}
