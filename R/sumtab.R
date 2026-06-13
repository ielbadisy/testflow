#' Build a compact descriptive summary table
#'
#' @param formula A one-sided formula such as `~ age + sex` or
#'   `~ age + sex | treatment`.
#' @param data A data frame.
#' @param p_value Logical; add a p-value column when a grouping variable is
#'   supplied.
#' @param overall Logical; include an overall summary column.
#' @param digits Number of digits for summary statistics.
#' @param p_digits Number of digits for formatted p-values.
#' @param alpha Significance level used by automatic test selection.
#' @param fisher_threshold Expected-count threshold for Fisher's exact test.
#' @param na.rm Logical; remove missing values before summaries and tests.
#' @return A tibble with one row per numeric variable and one row per
#'   categorical level.
#' @export
sumtab <- function(
  formula,
  data,
  p_value = FALSE,
  overall = TRUE,
  digits = 1,
  p_digits = 3,
  alpha = 0.05,
  fisher_threshold = 5,
  na.rm = TRUE
) {
  parsed <- parse_sumtab_formula(substitute(formula))
  require_columns(data, c(parsed$vars, parsed$group))

  group_levels <- sumtab_group_levels(data, parsed$group)

  purrr::map_dfr(parsed$vars, function(variable) {
    rows <- if (is.numeric(data[[variable]])) {
      sumtab_numeric_rows(data, variable, parsed$group, group_levels, overall, digits, na.rm)
    } else {
      sumtab_categorical_rows(data, variable, parsed$group, group_levels, overall, digits, na.rm)
    }

    if (!is.null(parsed$group) && isTRUE(p_value)) {
      test <- sumtab_auto_test(data, variable, parsed$group, alpha, fisher_threshold, na.rm)
      rows$p.value <- c(format_p(test$p.value, p_digits), rep(NA_character_, nrow(rows) - 1))
      rows$test <- c(test$method, rep(NA_character_, nrow(rows) - 1))
    }

    rows
  })
}

parse_sumtab_formula <- function(expr) {
  is_formula <- inherits(expr, "formula") || (is.call(expr) && identical(expr[[1]], as.name("~")))
  if (!is_formula || length(expr) != 2) {
    stop("`formula` must be one-sided, for example `~ age + sex | treatment`.", call. = FALSE)
  }

  rhs <- expr[[2]]
  group <- NULL
  vars_expr <- rhs

  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))) {
    vars_expr <- rhs[[2]]
    group <- all.vars(rhs[[3]])
    if (length(group) != 1) {
      stop("The grouping side of `formula` must contain exactly one variable.", call. = FALSE)
    }
  }

  vars <- all.vars(vars_expr)
  if (length(vars) == 0) {
    stop("`formula` must select at least one summary variable.", call. = FALSE)
  }

  list(vars = vars, group = group)
}

sumtab_group_levels <- function(data, group) {
  if (is.null(group)) {
    return(character())
  }

  x <- data[[group]]
  if (is.factor(x)) {
    return(as.character(levels(droplevels(x[!is.na(x)]))))
  }

  as.character(unique(stats::na.omit(x)))
}

sumtab_numeric_rows <- function(data, variable, group, group_levels, overall, digits, na.rm) {
  cells <- list()
  if (isTRUE(overall)) {
    cells[["Overall"]] <- sumtab_numeric_cell(data[[variable]], digits, na.rm)
  }
  for (level in group_levels) {
    cells[[level]] <- sumtab_numeric_cell(data[[variable]][data[[group]] == level], digits, na.rm)
  }

  tibble::tibble(variable = variable, level = NA_character_) |>
    dplyr::bind_cols(tibble::as_tibble_row(cells))
}

sumtab_categorical_rows <- function(data, variable, group, group_levels, overall, digits, na.rm) {
  levels <- sumtab_variable_levels(data[[variable]], na.rm)

  purrr::map_dfr(levels, function(level) {
    cells <- list()
    if (isTRUE(overall)) {
      cells[["Overall"]] <- sumtab_categorical_cell(data[[variable]], level, digits, na.rm)
    }
    for (group_level in group_levels) {
      cells[[group_level]] <- sumtab_categorical_cell(data[[variable]][data[[group]] == group_level], level, digits, na.rm)
    }

    tibble::tibble(variable = variable, level = level) |>
      dplyr::bind_cols(tibble::as_tibble_row(cells))
  })
}

sumtab_variable_levels <- function(x, na.rm) {
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  }
  if (is.factor(x)) {
    return(as.character(levels(droplevels(x))))
  }
  as.character(unique(x))
}

sumtab_numeric_cell <- function(x, digits, na.rm) {
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  }
  n <- length(x)
  if (n == 0) {
    return("0")
  }

  paste0(
    format_stat(mean(x), digits),
    " (",
    format_stat(stats::sd(x), digits),
    "); ",
    format_stat(stats::median(x), digits),
    " [",
    format_stat(stats::quantile(x, 0.25, names = FALSE), digits),
    ", ",
    format_stat(stats::quantile(x, 0.75, names = FALSE), digits),
    "]; n=",
    n
  )
}

sumtab_categorical_cell <- function(x, level, digits, na.rm) {
  if (isTRUE(na.rm)) {
    x <- x[!is.na(x)]
  }
  n_total <- length(x)
  n <- sum(as.character(x) == level, na.rm = TRUE)
  pct <- if (n_total == 0) NA_real_ else n / n_total * 100
  paste0(n, " (", format_stat(pct, digits), "%)")
}
