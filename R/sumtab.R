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
  group_labels <- sumtab_group_labels(data, parsed$group, group_levels)
  overall_label <- paste0("Overall (n = ", nrow(data), ")")

  purrr::map_dfr(parsed$vars, function(variable) {
    rows <- if (is.numeric(data[[variable]])) {
      sumtab_numeric_rows(data, variable, parsed$group, group_levels, group_labels, overall, overall_label, digits, na.rm)
    } else {
      sumtab_categorical_rows(data, variable, parsed$group, group_levels, group_labels, overall, overall_label, digits, na.rm)
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

sumtab_group_labels <- function(data, group, group_levels) {
  if (is.null(group)) {
    return(stats::setNames(character(), character()))
  }

  counts <- stats::setNames(
    vapply(group_levels, function(level) sum(data[[group]] == level, na.rm = TRUE), integer(1)),
    group_levels
  )
  stats::setNames(paste0(names(counts), " (n = ", counts, ")"), names(counts))
}

sumtab_numeric_rows <- function(data, variable, group, group_levels, group_labels, overall, overall_label, digits, na.rm) {
  cells <- list()
  if (isTRUE(overall)) {
    cells[[overall_label]] <- sumtab_numeric_cell(data[[variable]], digits, na.rm)
  }
  for (level in group_levels) {
    cells[[group_labels[[level]]]] <- sumtab_numeric_cell(data[[variable]][data[[group]] == level], digits, na.rm)
  }

  tibble::tibble(variable = variable, level = NA_character_) |>
    dplyr::bind_cols(tibble::as_tibble_row(cells))
}

sumtab_categorical_rows <- function(data, variable, group, group_levels, group_labels, overall, overall_label, digits, na.rm) {
  levels <- sumtab_variable_levels(data[[variable]], na.rm)

  purrr::map_dfr(levels, function(level) {
    cells <- list()
    if (isTRUE(overall)) {
      cells[[overall_label]] <- sumtab_categorical_cell(data[[variable]], level, digits, na.rm)
    }
    for (group_level in group_levels) {
      cells[[group_labels[[group_level]]]] <- sumtab_categorical_cell(data[[variable]][data[[group]] == group_level], level, digits, na.rm)
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

sumtab_auto_test <- function(data, variable, group, alpha, fisher_threshold, na.rm) {
  df <- drop_missing(data, c(variable, group), na.rm = na.rm)
  df[[group]] <- as.factor(df[[group]])
  n_groups <- dplyr::n_distinct(df[[group]])

  if (n_groups < 2) {
    return(list(p.value = NA_real_, method = NA_character_))
  }

  if (is.numeric(df[[variable]])) {
    if (n_groups == 2) {
      return(sumtab_numeric_two_group_test(df, variable, group, alpha))
    }
    return(sumtab_numeric_multi_group_test(df, variable, group, alpha))
  }

  sumtab_categorical_test(df, variable, group, fisher_threshold)
}

sumtab_numeric_two_group_test <- function(data, variable, group, alpha) {
  normality <- check_normality(data, variable, group, alpha)
  variance <- check_variance_homogeneity(data, variable, group, alpha)
  recommendation <- recommend_two_groups(normality, variance)
  formula <- stats::as.formula(paste(variable, "~", group))

  test <- switch(
    recommendation,
    "Student independent t-test" = stats::t.test(formula, data = data, var.equal = TRUE),
    "Welch t-test" = stats::t.test(formula, data = data, var.equal = FALSE),
    "Wilcoxon rank-sum test" = stats::wilcox.test(formula, data = data, exact = FALSE)
  )

  list(p.value = test$p.value, method = recommendation)
}

sumtab_numeric_multi_group_test <- function(data, variable, group, alpha) {
  normality <- check_normality(data, variable, group, alpha)
  variance <- check_variance_homogeneity(data, variable, group, alpha)
  recommendation <- recommend_groups(normality, variance)
  formula <- stats::as.formula(paste(variable, "~", group))

  p <- switch(
    recommendation,
    "One-way ANOVA" = stats::anova(stats::aov(formula, data = data))[["Pr(>F)"]][1],
    "Welch ANOVA" = stats::oneway.test(formula, data = data, var.equal = FALSE)$p.value,
    "Kruskal-Wallis test" = stats::kruskal.test(formula, data = data)$p.value
  )

  list(p.value = p, method = recommendation)
}

sumtab_categorical_test <- function(data, variable, group, fisher_threshold) {
  tab <- table(data[[variable]], data[[group]])
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))

  if (any(chi$expected < fisher_threshold)) {
    fisher <- stats::fisher.test(tab)
    return(list(p.value = fisher$p.value, method = "Fisher exact test"))
  }

  list(p.value = chi$p.value, method = "Chi-square test of independence")
}
