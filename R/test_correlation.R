#' Test correlation between two numeric variables
#' @param data A data frame.
#' @param x First numeric column.
#' @param y Second numeric column.
#' @param method Correlation method or `"auto"`.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @export
test_correlation <- function(data, x, y, method = c("auto", "pearson", "spearman", "kendall"), alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  method <- match.arg(method)
  x_nm <- rlang::as_name(rlang::ensym(x)); y_nm <- rlang::as_name(rlang::ensym(y))
  df <- drop_missing(data, c(x_nm, y_nm), na.rm = na.rm)
  normality <- check_normality(df, c(x_nm, y_nm), alpha = alpha)
  outliers <- iqr_outliers(df, c(x_nm, y_nm))
  chosen <- if (method == "auto") {
    if (all(normality$status == "acceptable") && !any(outliers$is_outlier)) "pearson" else "spearman"
  } else method
  pearson <- stats::cor.test(df[[x_nm]], df[[y_nm]], method = "pearson")
  spearman <- suppressWarnings(stats::cor.test(df[[x_nm]], df[[y_nm]], method = "spearman", exact = FALSE))
  kendall <- suppressWarnings(stats::cor.test(df[[x_nm]], df[[y_nm]], method = "kendall", exact = FALSE))
  primary <- switch(chosen, pearson = pearson, spearman = spearman, kendall = kendall)
  effect <- tibble::tibble(name = paste0(chosen, " r"), estimate = unname(primary$estimate), magnitude = magnitude_cramers_v(abs(unname(primary$estimate))))
  plt <- if (plot) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_nm]], y = .data[[y_nm]])) +
      ggplot2::geom_point(alpha = 0.75, color = "#4C78A8") +
      ggplot2::geom_smooth(method = "lm", se = TRUE, color = "#F58518") +
      ggplot2::labs(title = "Correlation workflow", subtitle = plot_subtitle(paste(chosen, "correlation"), primary), caption = paste0(effect$name[1], " = ", format_stat(effect$estimate[1]), ", ", effect$magnitude[1]), x = x_nm, y = y_nm) +
      ggplot2::theme_minimal()
  } else NULL
  out <- new_testflow("correlation", "two numeric variables", y_nm, x_nm, data = df, descriptives = descriptives_numeric(df, c(x_nm, y_nm)), assumptions = list("Normality" = normality, "IQR outliers" = outliers), recommended = list(test = paste(chosen, "correlation")), primary_test = safe_tidy_htest(primary, paste(chosen, "correlation")), alternative_tests = list(pearson = safe_tidy_htest(pearson, "Pearson correlation"), spearman = safe_tidy_htest(spearman, "Spearman correlation"), kendall = safe_tidy_htest(kendall, "Kendall correlation")), effect_size = effect, plot = plt, call = match.call(), subclass = "correlation")
  out$interpretation <- make_report(out, alpha)
  out
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
  df <- drop_missing(data, vars, na.rm = na.rm)
  mat <- stats::cor(df[, vars, drop = FALSE], method = method, use = "pairwise.complete.obs")
  pairs <- utils::combn(vars, 2, simplify = FALSE)
  long <- purrr::map_dfr(pairs, function(z) {
    test <- suppressWarnings(stats::cor.test(df[[z[1]]], df[[z[2]]], method = method, exact = FALSE))
    tibble::tibble(var1 = z[1], var2 = z[2], estimate = unname(test$estimate), p = test$p.value)
  })
  plt <- if (plot) {
    hm <- as.data.frame(as.table(mat))
    ggplot2::ggplot(hm, ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$Freq)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
      ggplot2::labs(title = "Correlation matrix workflow", x = NULL, y = NULL, fill = "r") +
      ggplot2::theme_minimal()
  } else NULL
  primary <- tibble::tibble(method = paste(method, "correlation matrix"), statistic = NA_real_, parameter = NA_real_, p.value = min(long$p, na.rm = TRUE))
  out <- new_testflow("correlation_matrix", "correlation matrix", paste(vars, collapse = ", "), data = df, descriptives = descriptives_numeric(df, vars), recommended = list(test = paste(method, "correlation matrix")), primary_test = primary, alternative_tests = list(correlation_matrix = mat, p_values = long), effect_size = tibble::tibble(name = "maximum absolute r", estimate = max(abs(long$estimate), na.rm = TRUE), magnitude = magnitude_cramers_v(max(abs(long$estimate), na.rm = TRUE))), plot = plt, call = match.call(), subclass = "correlation_matrix")
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
