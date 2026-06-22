#' Run a repeated-measures workflow from wide data
#' @param data A data frame.
#' @param measures Repeated numeric columns selected with tidyselect syntax.
#' @param id Optional subject identifier.
#' @param between Optional between-subject factor.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @references
#' Fisher, R. A. (1925). \emph{Statistical Methods for Research Workers}.
#' Oliver and Boyd.
#'
#' Friedman, M. (1937). The use of ranks to avoid the assumption of normality
#' implicit in the analysis of variance. \emph{Journal of the American
#' Statistical Association}, 32(200), 675-701.
#'
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods.
#' \emph{Biometrics Bulletin}, 1(6), 80-83.
#'
#' Girden, E. R. (1992). \emph{ANOVA: Repeated Measures}. Sage.
#' @export
test_repeated <- function(data, measures, id = NULL, between = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  measure_nms <- tidyselect_names(data, {{ measures }})
  warn_if_two_measures_repeated(measure_nms)
  id_nm <- if (missing(id) || is.null(substitute(id))) ".testflow_id" else rlang::as_name(rlang::ensym(id))
  df <- data
  if (id_nm == ".testflow_id") df[[id_nm]] <- seq_len(nrow(df))
  long <- tidyr::pivot_longer(df, dplyr::all_of(measure_nms), names_to = "time", values_to = "value")
  repeated_core(long, "value", "time", id_nm, outcome_label = paste(measure_nms, collapse = ", "), alpha = alpha, plot = plot, na.rm = na.rm, call = match.call())
}

#' Run a repeated-measures workflow from long data
#' @param data A data frame.
#' @param outcome Numeric outcome column.
#' @param within Within-subject time/condition column.
#' @param id Subject identifier column.
#' @param between Optional between-subject factor.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_repeated_long <- function(data, outcome, within, id, between = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  within_nm <- rlang::as_name(rlang::ensym(within))
  id_nm <- rlang::as_name(rlang::ensym(id))
  repeated_core(data, outcome_nm, within_nm, id_nm, outcome_label = outcome_nm, alpha = alpha, plot = plot, na.rm = na.rm, call = match.call())
}

repeated_core <- function(data, outcome_nm, within_nm, id_nm, outcome_label = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE, call = NULL) {
  if (is.null(outcome_label)) outcome_label <- outcome_nm
  df <- drop_missing(data, c(outcome_nm, within_nm, id_nm), na.rm = na.rm)
  normality <- check_normality(df, outcome_nm, within_nm, alpha)
  sphericity <- check_sphericity_or_note()
  wide <- tidyr::pivot_wider(df, names_from = dplyr::all_of(within_nm), values_from = dplyr::all_of(outcome_nm), id_cols = dplyr::all_of(id_nm))
  measure_nms <- setdiff(names(wide), id_nm)
  complete <- stats::na.omit(wide[, measure_nms, drop = FALSE])
  friedman <- stats::friedman.test(as.matrix(complete))
  anova <- repeated_anova_test(df, outcome_nm, within_nm, id_nm)
  recommendation <- if (all(normality$status == "acceptable")) "Repeated-measures ANOVA" else "Friedman test"
  h0 <- h0_mean_equal(outcome_nm, within_nm)
  primary <- if (recommendation == "Friedman test") {
    add_null_hypothesis(safe_tidy_htest(friedman, "Friedman test"), h0)
  } else {
    add_null_hypothesis(anova, h0)
  }
  effect <- if (recommendation == "Friedman test") {
    friedman_effect_size(friedman, n = nrow(complete), k = length(measure_nms))
  } else {
    repeated_eta_squared(df, outcome_nm, within_nm, id_nm)
  }
  posthoc_method <- if (recommendation == "Friedman test") "wilcox" else "t"
  posthoc_tests <- paired_posthoc_numeric(complete, method = posthoc_method, p.adjust.method = "BH")
  plt <- if (plot) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[within_nm]], y = .data[[outcome_nm]], group = .data[[id_nm]])) +
      ggplot2::geom_line(alpha = 0.25) +
      ggplot2::geom_point(alpha = 0.45) +
      ggplot2::stat_summary(ggplot2::aes(group = 1), fun = mean, geom = "line", linewidth = 1.1, color = "#F58518") +
      ggplot2::labs(title = "Repeated-measures workflow", subtitle = paste0(recommendation, ", p = ", format_p(primary$p.value[1])), x = within_nm, y = outcome_nm) +
      ggplot2::theme_minimal()
  } else NULL
  out <- new_testflow("repeated", "repeated numeric measurements", outcome_label, within_nm, id_nm, df, descriptives_numeric(df, outcome_nm, within_nm), assumption_checks(check_independence_note("Repeated measurements from the same subjects are assumed by design."), normality, sphericity), list(test = recommendation), primary, list(anova = add_null_hypothesis(anova, h0), friedman = add_null_hypothesis(safe_tidy_htest(friedman, "Friedman test"), h0)), posthoc = posthoc_tests, effect_size = effect, plot = plt, call = call, subclass = "repeated")
  out$interpretation <- make_report(out, alpha)
  out
}
