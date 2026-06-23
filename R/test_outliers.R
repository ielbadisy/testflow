#' Detect numeric outliers
#' @param formula Numeric columns to screen for outliers.
#' @param data A data frame.
#' @param group Optional grouping column.
#' @param method Outlier method.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_outliers`. The object is a
#' list containing the cleaned data, numeric descriptives, screening
#' assumptions, selected outlier-screening method, flagged IQR and/or
#' Mahalanobis rows, outlier-count summary, optional `ggplot`, original call,
#' and report text. This is a screening workflow, not a single hypothesis test.
#' @export
test_outliers <- function(formula, data, group = NULL, method = c("iqr", "mahalanobis", "both"), plot = TRUE, na.rm = TRUE) {
  method <- match.arg(method)
  vars <- tidyselect_names(data, {{ formula }})
  group_nm <- if (missing(group) || is.null(substitute(group))) NULL else rlang::as_name(rlang::ensym(group))
  warn_if_screening_workflow("outliers")
  warn_if(method == "mahalanobis" && length(vars) < 2, "Mahalanobis outlier screening requires at least two numeric variables.")
  cols <- c(vars, group_nm)
  df <- drop_missing(data, cols, na.rm = na.rm)
  iqr <- if (method %in% c("iqr", "both")) iqr_outliers(df, vars, group_nm) else NULL
  mahal <- if (method %in% c("mahalanobis", "both") && length(vars) > 1) mahalanobis_outliers(df, vars) else NULL
  assumptions <- if (method == "iqr") {
    assumption_checks(assumption_check("Numeric variable", "acceptable", "IQR outlier detection is univariate and does not require normality."), assumption_check("Skewness sensitivity", "warning", "Interpret IQR outliers with care when the distribution is strongly skewed."))
  } else {
    assumption_checks(assumption_check("Complete cases", "acceptable", "Mahalanobis distance uses complete cases only."), assumption_check("Approximate multivariate normality", "warning", "Mahalanobis screening is more stable when the variables are roughly multivariate normal."), assumption_check("Invertible covariance matrix", "acceptable", "The covariance matrix must be invertible."))
  }
  method_label <- ifelse(method == "iqr", "IQR outlier detection", title_case_method(paste(method, "outlier detection")))
  primary <- tibble::tibble(method = method_label, statistic = sum((iqr$is_outlier %||% FALSE), na.rm = TRUE), parameter = NA_real_, p.value = NA_real_)
  plt <- if (plot) {
    long <- tidyr::pivot_longer(df, dplyr::all_of(vars), names_to = "variable", values_to = "value")
    if (!is.null(iqr)) {
      iqr_flags <- dplyr::select(iqr, "row", "variable", "is_outlier", "is_extreme")
      long <- dplyr::left_join(
        dplyr::mutate(long, row = dplyr::row_number(), .by = "variable"),
        iqr_flags,
        by = c("row", "variable")
      )
    }
    fences <- iqr_fences(df, vars)
    ggplot2::ggplot(long, ggplot2::aes(x = .data$variable, y = .data$value)) +
      ggplot2::geom_hline(data = fences, ggplot2::aes(yintercept = .data$lower_fence), linetype = "dashed", color = "#B23A48") +
      ggplot2::geom_hline(data = fences, ggplot2::aes(yintercept = .data$upper_fence), linetype = "dashed", color = "#B23A48") +
      ggplot2::geom_boxplot(fill = "#D9E2EC", color = "#334E68", width = 0.38, outlier.shape = NA) +
      ggplot2::geom_jitter(ggplot2::aes(color = .data$is_outlier), width = 0.08, alpha = 0.7, size = 1.6, na.rm = TRUE) +
      ggplot2::scale_color_manual(values = c(`FALSE` = "#52606D", `TRUE` = "#D64545"), breaks = c(FALSE, TRUE), labels = c("inside fence", "IQR outlier"), na.translate = FALSE) +
      ggplot2::facet_wrap(~variable, scales = "free_y", nrow = 1) +
      ggplot2::labs(title = "Outlier screening", subtitle = "Dashed lines mark Q1 - 1.5 x IQR and Q3 + 1.5 x IQR", x = NULL, y = NULL, color = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank())
  } else NULL
  out <- new_testflow("outliers", "outlier screening", paste(vars, collapse = ", "), group_nm, data = df, descriptives = descriptives_numeric(df, vars, group_nm), assumptions = assumptions, recommended = list(test = method_label), primary_test = primary, alternative_tests = list(iqr = iqr, mahalanobis = mahal), effect_size = tibble::tibble(name = "Outlier count", estimate = primary$statistic[1], magnitude = NA_character_), plot = plt, call = match.call(), subclass = "outliers")
  out$interpretation <- make_report(out)
  out
}

iqr_fences <- function(data, vars) {
  purrr::map_dfr(vars, function(v) {
    q1 <- stats::quantile(data[[v]], 0.25, na.rm = TRUE, names = FALSE)
    q3 <- stats::quantile(data[[v]], 0.75, na.rm = TRUE, names = FALSE)
    i <- q3 - q1
    tibble::tibble(variable = v, lower_fence = q1 - 1.5 * i, upper_fence = q3 + 1.5 * i)
  })
}

iqr_outliers <- function(data, vars, group = NULL) {
  calc <- function(df) {
    purrr::map_dfr(vars, function(v) {
      q1 <- stats::quantile(df[[v]], 0.25, na.rm = TRUE)
      q3 <- stats::quantile(df[[v]], 0.75, na.rm = TRUE)
      i <- q3 - q1
      tibble::tibble(row = seq_len(nrow(df)), variable = v, value = df[[v]], is_outlier = df[[v]] < q1 - 1.5 * i | df[[v]] > q3 + 1.5 * i, is_extreme = df[[v]] < q1 - 3 * i | df[[v]] > q3 + 3 * i)
    })
  }
  if (is.null(group)) calc(data) else data |> dplyr::group_by(.data[[group]]) |> dplyr::group_modify(~ calc(.x)) |> dplyr::ungroup()
}

mahalanobis_outliers <- function(data, vars) {
  x <- stats::na.omit(as.matrix(data[, vars, drop = FALSE]))
  d <- stats::mahalanobis(x, colMeans(x), stats::cov(x))
  cutoff <- stats::qchisq(0.975, df = length(vars))
  tibble::tibble(row = seq_along(d), distance = d, cutoff = cutoff, is_outlier = d > cutoff)
}
