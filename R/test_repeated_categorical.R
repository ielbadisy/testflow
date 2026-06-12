#' Test repeated categorical measurements
#' @param data A data frame.
#' @param measures Repeated binary columns selected with tidyselect syntax.
#' @param id Optional subject identifier.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_repeated_categorical <- function(data, measures, id = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  measure_nms <- tidyselect_names(data, {{ measures }})
  df <- drop_missing(data, measure_nms, na.rm = na.rm)
  mat <- as.matrix(df[, measure_nms, drop = FALSE])
  mat <- apply(mat, 2, function(x) as.integer(as.factor(x)) - 1L)
  test <- cochran_q_test(mat)
  counts <- tibble::as_tibble(mat) |> dplyr::summarise(dplyr::across(dplyr::everything(), sum)) |> tidyr::pivot_longer(dplyr::everything(), names_to = "time", values_to = "success")
  effect <- cochran_q_effect_size(test, n = nrow(mat), k = ncol(mat))
  plt <- if (plot) {
    counts$prop <- counts$success / nrow(mat)
    ggplot2::ggplot(counts, ggplot2::aes(x = .data$time, y = .data$prop, group = 1)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(title = "Repeated categorical workflow", subtitle = plot_subtitle("Cochran Q test", test), x = NULL, y = "Proportion") +
      ggplot2::theme_minimal()
  } else NULL
  h0 <- paste0("H0: the success proportions are equal across repeated categorical measures.")
  out <- new_testflow("repeated_categorical", "repeated categorical measurements", paste(measure_nms, collapse = ", "), data = df, descriptives = counts, recommended = list(test = "Cochran Q test"), primary_test = add_null_hypothesis(safe_tidy_htest(test, "Cochran Q test"), h0), alternative_tests = list(pairwise_mcnemar = pairwise_mcnemar(mat, measure_nms)), posthoc = pairwise_mcnemar(mat, measure_nms), effect_size = effect, plot = plt, call = match.call(), subclass = "repeated_categorical")
  out$interpretation <- make_report(out, alpha)
  out
}
