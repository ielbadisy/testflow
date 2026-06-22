#' Test correlation between two numeric variables
#' @param formula A formula such as `y ~ x`, or a data frame when using pipe/data-first style.
#' @param data A data frame, or a first numeric column when using data-first style.
#' @param y Second numeric column. Optional when using formula style.
#' @param method Correlation method or `"auto"`.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @references
#' Pearson, K. (1895). Notes on regression and inheritance in the case of two
#' parents. \emph{Proceedings of the Royal Society of London}, 58, 240-242.
#'
#' Spearman, C. (1904). The proof and measurement of association between two
#' things. \emph{The American Journal of Psychology}, 15(1), 72-101.
#'
#' Kendall, M. G. (1938). A new measure of rank correlation. \emph{Biometrika},
#' 30(1/2), 81-93.
#' @export
test_correlation <- function(formula, data, y = NULL, method = c("auto", "pearson", "spearman", "kendall"), alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  method <- match.arg(method)
  vars <- resolve_formula_pair(substitute(formula), substitute(data), substitute(y), missing(y), rhs_label = "x")
  y_nm <- vars$outcome
  x_nm <- vars$group
  data_obj <- resolve_data_first_or_formula(formula, data)
  df <- drop_missing(data_obj, c(x_nm, y_nm), na.rm = na.rm)
  warn_if_small_correlation_n(df, x_nm, y_nm)
  normality <- check_normality(df, c(x_nm, y_nm), alpha = alpha)
  outlier_flags <- iqr_outlier_flags(df, c(x_nm, y_nm))
  outliers <- iqr_outlier_assumption(outlier_flags, c(x_nm, y_nm))
  chosen <- if (method == "auto") {
    if (all(normality$status == "acceptable") && !any(outlier_flags$is_outlier, na.rm = TRUE)) "pearson" else "spearman"
  } else method
  linearity <- assumption_check("Linearity", ifelse(chosen == "pearson", "acceptable", "not checked"), ifelse(chosen == "pearson", "A roughly linear relation is assumed for Pearson correlation.", "Normality is not required for Spearman or Kendall; monotonicity is the key check."))
  monotonicity <- check_monotonicity(df[[x_nm]], df[[y_nm]], alpha = alpha)
  pearson <- stats::cor.test(df[[x_nm]], df[[y_nm]], method = "pearson")
  spearman <- suppressWarnings(stats::cor.test(df[[x_nm]], df[[y_nm]], method = "spearman", exact = FALSE))
  kendall <- suppressWarnings(stats::cor.test(df[[x_nm]], df[[y_nm]], method = "kendall", exact = FALSE))
  primary <- switch(chosen, pearson = pearson, spearman = spearman, kendall = kendall)
  chosen_label <- title_case_method(paste(chosen, "correlation"))
  effect <- tibble::tibble(name = paste0(chosen_label, " r"), estimate = unname(primary$estimate), magnitude = magnitude_cramers_v(abs(unname(primary$estimate))))
  plt <- if (plot) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_nm]], y = .data[[y_nm]])) +
      ggplot2::geom_point(alpha = 0.75, color = "#4C78A8") +
      ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#F58518") +
      ggplot2::labs(title = "Correlation workflow", subtitle = plot_subtitle(chosen_label, primary), caption = paste0(effect$name[1], " = ", format_stat(effect$estimate[1]), ", ", effect$magnitude[1]), x = x_nm, y = y_nm) +
      ggplot2::theme_minimal()
  } else NULL
  h0 <- h0_no_correlation(x_nm, y_nm)
  corr_assumptions <- switch(chosen, pearson = assumption_checks(linearity, normality, outliers), spearman = assumption_checks(monotonicity, outliers, assumption_check("Normality", "not required", "Normality is not required for Spearman correlation.")), kendall = assumption_checks(monotonicity, outliers, assumption_check("Normality", "not required", "Normality is not required for Kendall correlation.")))
  out <- new_testflow("correlation", "two numeric variables", y_nm, x_nm, data = df, descriptives = descriptives_numeric(df, c(x_nm, y_nm)), assumptions = corr_assumptions, recommended = list(test = chosen_label), primary_test = add_null_hypothesis(safe_tidy_htest(primary, chosen_label), h0), alternative_tests = list(correlation_table = tibble::tibble(method = c("Pearson", "Spearman", "Kendall"), statistic = c(unname(pearson$estimate), unname(spearman$estimate), unname(kendall$estimate)), p.value = c(pearson$p.value, spearman$p.value, kendall$p.value)), pearson = add_null_hypothesis(safe_tidy_htest(pearson, "Pearson correlation"), h0), spearman = add_null_hypothesis(safe_tidy_htest(spearman, "Spearman correlation"), h0), kendall = add_null_hypothesis(safe_tidy_htest(kendall, "Kendall correlation"), h0)), effect_size = effect, plot = plt, call = match.call(), subclass = "correlation")
  out$interpretation <- make_report(out, alpha)
  out
}

iqr_outlier_assumption <- function(data, vars) {
  n <- sum(data$is_outlier, na.rm = TRUE)
  assumption_check(
    "Extreme outliers",
    ifelse(n == 0, "acceptable", "warning"),
    paste0(n, " potential outlier(s) flagged by IQR."),
    details = paste0("IQR rule applied to ", paste(vars, collapse = ", "))
  )
}

iqr_outlier_flags <- function(data, vars) {
  purrr::map_dfr(vars, function(v) {
    x <- data[[v]]
    x <- x[!is.na(x)]
    if (length(x) < 4) {
      return(tibble::tibble(variable = v, is_outlier = NA))
    }
    q1 <- stats::quantile(x, 0.25, names = FALSE)
    q3 <- stats::quantile(x, 0.75, names = FALSE)
    i <- q3 - q1
    tibble::tibble(variable = v, is_outlier = x < q1 - 1.5 * i | x > q3 + 1.5 * i)
  })
}

#' Test a correlation matrix
#' @param data A data frame.
#' @param vars Numeric columns.
#' @param method Correlation method.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_correlation_matrix <- function(data, vars, method = c("spearman", "pearson", "kendall"), alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  method <- match.arg(method)
  vars <- tidyselect_names(data, {{ vars }})
  warn_if_screening_workflow("correlation_matrix")
  df <- drop_missing(data, vars, na.rm = na.rm)
  mat <- stats::cor(df[, vars, drop = FALSE], method = method, use = "pairwise.complete.obs")
  pairs <- utils::combn(vars, 2, simplify = FALSE)
  long <- purrr::map_dfr(pairs, function(z) {
    test <- suppressWarnings(stats::cor.test(df[[z[1]]], df[[z[2]]], method = method, exact = FALSE))
    tibble::tibble(var1 = z[1], var2 = z[2], estimate = unname(test$estimate), p = test$p.value)
  })
  assumptions <- assumption_checks(
    assumption_check("Variable types", "acceptable", "All selected variables should be numeric or ordinal."),
    assumption_check("Pairwise missing data", "acceptable", "Pairwise complete observations are used for each correlation."),
    assumption_check("Multiple testing correction", "not checked", "P-values are reported pairwise; correction is not applied in this workflow.")
  )
  plt <- if (plot) {
    hm <- as.data.frame(as.table(mat))
    ggplot2::ggplot(hm, ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$Freq)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
      ggplot2::labs(title = "Correlation matrix workflow", x = NULL, y = NULL, fill = "r") +
      ggplot2::theme_minimal()
  } else NULL
  method_label <- title_case_method(paste(method, "correlation matrix"))
  primary <- tibble::tibble(method = method_label, statistic = NA_real_, parameter = NA_real_, p.value = min(long$p, na.rm = TRUE))
  out <- new_testflow("correlation_matrix", "correlation matrix", paste(vars, collapse = ", "), data = df, descriptives = descriptives_numeric(df, vars), assumptions = assumptions, recommended = list(test = method_label), primary_test = primary, alternative_tests = list(correlation_matrix = mat, correlation_table = long), effect_size = tibble::tibble(name = "Maximum absolute r", estimate = max(abs(long$estimate), na.rm = TRUE), magnitude = magnitude_cramers_v(max(abs(long$estimate), na.rm = TRUE))), plot = plt, call = match.call(), subclass = "correlation_matrix")
  out$interpretation <- make_report(out, alpha)
  out
}

tidyselect_names <- function(data, expr) {
  if (rlang::quo_is_symbol(rlang::enquo(expr))) {
    nm <- rlang::as_name(rlang::ensym(expr))
    if (nm %in% names(data)) return(nm)
  }
  tidyselect::eval_select(rlang::enquo(expr), data) |> names()
}
