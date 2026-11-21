#' Test paired categorical measurements
#' @param data A data frame.
#' @param before Before categorical column.
#' @param after After categorical column.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param correct Logical; apply the continuity correction in McNemar's test,
#'   \eqn{\chi^2 = (|b - c| - 1)^2 / (b + c)}, where `b` and `c` are the
#'   discordant counts. Defaults to `TRUE` (as in [stats::mcnemar.test()]);
#'   `FALSE` gives the uncorrected \eqn{(b - c)^2 / (b + c)}.
#' @return A `testflow` object with class `testflow_paired_categorical`. The
#' object is a list containing the cleaned paired categorical data, categorical
#' descriptives, assumption checks, McNemar test result, discordant-pair table,
#' optional `ggplot`, original call, and report text.
#' @export
testpairedcategorical <- function(data, before, after, alpha = 0.05, plot = TRUE, na.rm = TRUE, correct = TRUE) {
  before_nm <- rlang::as_name(rlang::ensym(before)); after_nm <- rlang::as_name(rlang::ensym(after))
  df <- drop_missing(data, c(before_nm, after_nm), na.rm = na.rm)
  tab <- table(df[[before_nm]], df[[after_nm]])
  warn_if(!all(dim(tab) == 2), "McNemar's test requires a 2x2 table; the current input is not a binary paired design.")
  test <- stats::mcnemar.test(tab, correct = correct)
  discordant <- if (all(dim(tab) == 2)) tab[1, 2] + tab[2, 1] else NA_integer_
  disc <- .bt(before_only = if (all(dim(tab) == 2)) tab[2, 1] else NA_integer_, after_only = if (all(dim(tab) == 2)) tab[1, 2] else NA_integer_)
  effect <- .bt(name = "Discordant pairs", estimate = sum(disc, na.rm = TRUE), magnitude = NA_character_)
  plt <- if (plot) {
    # The 2 x 2 table of pairs: McNemar's test uses only the two off-diagonal (discordant) cells.
    td <- as.data.frame(tab, stringsAsFactors = FALSE); names(td) <- c("b", "a", "n")
    td$b <- factor(td$b, rev(rownames(tab))); td$a <- factor(td$a, colnames(tab))
    td$cell <- ifelse(as.character(td$b) == as.character(td$a), "same", ifelse(as.integer(factor(td$b, rownames(tab))) <
                      as.integer(td$a), "up", "down"))
    ggplot2::ggplot(td, ggplot2::aes(x = .data$a, y = .data$b, fill = .data$cell)) +
      ggplot2::geom_tile(colour = "white", linewidth = 1) +
      ggplot2::geom_text(ggplot2::aes(label = .data$n, colour = .data$cell == "same"), size = 4.5) +
      ggplot2::scale_fill_manual(values = c(same = "grey92", up = TFCOL[["a1"]], down = TFCOL[["a2"]]), guide = "none") +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "grey15", `FALSE` = "white"), guide = "none") +
      ggplot2::coord_fixed() +
      ggplot2::labs(subtitle = plot_subtitle("McNemar test", test),
                    caption = "Coloured cells: discordant pairs (the only ones McNemar's test uses).",
                    x = after_nm, y = before_nm) +
      tftheme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
  } else NULL
  discordant_message <- if (is.na(discordant)) {
    "Discordant pairs could not be counted."
  } else if (discordant >= 10) {
    "Discordant pairs are adequate for the McNemar approximation."
  } else {
    "Small numbers of discordant pairs favor exact McNemar."
  }
  out <- new_testflow("paired_categorical", "paired categorical measurements", paste(before_nm, after_nm, sep = " -> "), data = df, descriptives = descriptives_categorical(df, c(before_nm, after_nm)), assumptions = assumption_checks(assumption_check("Paired binary measurements", "assumed", "Same subjects should be measured twice."), assumption_check("Discordant pairs", ifelse(is.na(discordant) || discordant >= 10, "acceptable", "warning"), discordant_message, details = paste0("discordant = ", discordant))), recommended = list(test = "McNemar test"), primary_test = safe_tidy_htest(test, "McNemar test"), alternative_tests = list(discordant_pairs = disc), effect_size = effect, plot = plt, call = match.call(), subclass = "paired_categorical")
  out$interpretation <- make_report(out, alpha)
  out
}
