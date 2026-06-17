#' Test association between two categorical variables
#' @param formula A formula such as `x ~ y`, or a data frame when using pipe/data-first style.
#' @param data A data frame, or a first categorical column when using data-first style.
#' @param y Second categorical column. Optional when using formula style.
#' @param alpha Significance level.
#' @param fisher_threshold Expected-count threshold for Fisher's exact test.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
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
#' @export
test_categorical <- function(formula, data, y = NULL, alpha = 0.05, fisher_threshold = 5, plot = TRUE, na.rm = TRUE) {
  vars <- resolve_formula_pair(substitute(formula), substitute(data), substitute(y), missing(y), rhs_label = "y")
  x_nm <- vars$outcome
  y_nm <- vars$group
  data_obj <- resolve_data_first_or_formula(formula, data)
  df <- drop_missing(data_obj, c(x_nm, y_nm), na.rm = na.rm)
  tab <- table(df[[x_nm]], df[[y_nm]])
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  fisher <- stats::fisher.test(tab)
  recommendation <- if (any(chi$expected < fisher_threshold)) "Fisher exact test" else "Chi-square test of independence"
  primary <- if (recommendation == "Fisher exact test") fisher else chi
  effect <- cramers_v(tab)
  expected <- check_expected_counts(tab, fisher_threshold)
  plt <- if (plot) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_nm]], fill = .data[[y_nm]])) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(title = "Categorical association workflow", subtitle = plot_subtitle(recommendation, primary), caption = paste0(effect$name[1], " = ", format_stat(effect$estimate[1]), ", ", effect$magnitude[1]), x = x_nm, y = "Proportion", fill = y_nm) +
      ggplot2::theme_minimal()
  } else NULL
  h0 <- h0_no_association(x_nm, y_nm)
  out <- new_testflow("categorical", "two categorical variables", x_nm, y_nm, data = df, descriptives = descriptives_categorical(df, c(x_nm, y_nm)), assumptions = assumption_checks(check_independence_note(), expected), recommended = list(test = recommendation), primary_test = add_null_hypothesis(safe_tidy_htest(primary, recommendation), h0), alternative_tests = list(chi_square = add_null_hypothesis(safe_tidy_htest(chi, "Chi-square test of independence"), h0), fisher = add_null_hypothesis(safe_tidy_htest(fisher, "Fisher exact test"), h0)), effect_size = effect, plot = plt, call = match.call(), subclass = "categorical")
  out$interpretation <- make_report(out, alpha)
  out
}
