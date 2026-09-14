#' Numeric descriptive statistics
#' @param data A data frame.
#' @param vars Numeric columns.
#' @param group Optional grouping column.
#' @noRd
.descriptives_numeric_one <- function(data, vars, group_label) {
  .map_dfr(vars, function(v) {
    x <- data[[v]]
    tibble::tibble(
      variable = v,
      group = group_label,
      n = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE),
      iqr = stats::IQR(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  })
}

descriptives_numeric <- function(data, vars, group = NULL) {
  vars <- as.character(vars)
  if (is.null(group)) {
    return(.descriptives_numeric_one(data, vars, NA_character_))
  }

  g <- data[[group]]
  # dplyr::group_by() treats NA as its own group; match that (NA group last,
  # matching dplyr's sort order).
  lvls <- sort(unique(as.character(g[!is.na(g)])))
  if (anyNA(g)) lvls <- c(lvls, NA_character_)

  out <- .map_dfr(lvls, function(lv) {
    subset_rows <- if (is.na(lv)) is.na(g) else !is.na(g) & as.character(g) == lv
    .descriptives_numeric_one(data[subset_rows, , drop = FALSE], vars, lv)
  })
  # matches the original dplyr pipeline's column order for the grouped path
  # (group first), which differs from the ungrouped path (variable first).
  out[, c("group", setdiff(names(out), "group"))]
}

#' Categorical descriptive statistics
#' @param data A data frame.
#' @param vars Categorical columns.
#' @noRd
descriptives_categorical <- function(data, vars) {
  vars <- as.character(vars)
  .map_dfr(vars, function(v) {
    tab <- table(data[[v]], useNA = "ifany")
    tibble::tibble(
      variable = v,
      level = names(tab),
      n = as.integer(tab),
      percent = as.numeric(tab) / sum(tab) * 100
    )
  })
}
