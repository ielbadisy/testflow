#' Numeric descriptive statistics
#' @param data A data frame.
#' @param vars Numeric columns.
#' @param group Optional grouping column.
#' @noRd
descriptives_numeric <- function(data, vars, group = NULL) {
  vars <- as.character(vars)
  if (is.null(group)) {
    purrr::map_dfr(vars, function(v) {
      x <- data[[v]]
      tibble::tibble(
        variable = v,
        group = NA_character_,
        n = sum(!is.na(x)),
        mean = mean(x, na.rm = TRUE),
        sd = stats::sd(x, na.rm = TRUE),
        median = stats::median(x, na.rm = TRUE),
        iqr = stats::IQR(x, na.rm = TRUE),
        min = min(x, na.rm = TRUE),
        max = max(x, na.rm = TRUE)
      )
    })
  } else {
    data |>
      dplyr::group_by(.data[[group]]) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(vars),
          list(
            n = ~ sum(!is.na(.x)),
            mean = ~ mean(.x, na.rm = TRUE),
            sd = ~ stats::sd(.x, na.rm = TRUE),
            median = ~ stats::median(.x, na.rm = TRUE),
            iqr = ~ stats::IQR(.x, na.rm = TRUE),
            min = ~ min(.x, na.rm = TRUE),
            max = ~ max(.x, na.rm = TRUE)
          ),
          .names = "{.col}__{.fn}"
        ),
        .groups = "drop"
      ) |>
      tidyr::pivot_longer(-dplyr::all_of(group), names_to = "name", values_to = "value") |>
      tidyr::separate("name", into = c("variable", "stat"), sep = "__") |>
      tidyr::pivot_wider(names_from = "stat", values_from = "value") |>
      dplyr::rename(group = !!rlang::sym(group))
  }
}

#' Categorical descriptive statistics
#' @param data A data frame.
#' @param vars Categorical columns.
#' @noRd
descriptives_categorical <- function(data, vars) {
  vars <- as.character(vars)
  purrr::map_dfr(vars, function(v) {
    tab <- table(data[[v]], useNA = "ifany")
    tibble::tibble(
      variable = v,
      level = names(tab),
      n = as.integer(tab),
      percent = as.numeric(tab) / sum(tab) * 100
    )
  })
}
