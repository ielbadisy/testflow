as_name_quo <- function(x) {
  expr <- substitute(x, parent.frame())
  rlang::as_name(expr)
}

require_columns <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(data)
}

drop_missing <- function(data, cols, na.rm = TRUE) {
  require_columns(data, cols)
  if (!na.rm) {
    return(data)
  }
  stats::na.omit(data[, cols, drop = FALSE])
}

assert_two_groups <- function(data, group) {
  groups <- unique(data[[group]])
  groups <- groups[!is.na(groups)]
  if (length(groups) != 2) {
    stop("`group` must contain exactly two non-missing groups.", call. = FALSE)
  }
  groups
}

safe_tidy_htest <- function(x, method = NULL) {
  out <- suppressMessages(broom::tidy(x))
  out <- standardize_result_columns(out)
  if (!is.null(method)) {
    out$method <- method
  }
  out
}

add_null_hypothesis <- function(x, h0) {
  x <- standardize_result_columns(x)
  x$null_hypothesis <- h0
  x
}

standardize_result_columns <- function(x) {
  if (is.null(x) || !inherits(x, "data.frame")) return(x)
  if (!"statistic" %in% names(x)) x$statistic <- NA_real_
  if (!"parameter" %in% names(x)) x$parameter <- NA_real_
  if (!"p.value" %in% names(x)) x$p.value <- NA_real_
  if (!"conf.low" %in% names(x)) x$conf.low <- NA_real_
  if (!"conf.high" %in% names(x)) x$conf.high <- NA_real_
  if (!"null_hypothesis" %in% names(x)) x$null_hypothesis <- NA_character_
  if (!"effect_size" %in% names(x)) x$effect_size <- NA_real_
  if (!"effect_size_name" %in% names(x)) x$effect_size_name <- NA_character_
  if (!"effect_size_magnitude" %in% names(x)) x$effect_size_magnitude <- NA_character_
  x
}

add_effect_size_result <- function(x, effect_size) {
  x <- standardize_result_columns(x)
  if (!is.null(effect_size) && inherits(effect_size, "data.frame") && nrow(effect_size) > 0) {
    x$effect_size <- effect_size$estimate[1]
    x$effect_size_name <- effect_size$name[1]
    x$effect_size_magnitude <- effect_size$magnitude[1]
  }
  x
}

h0_mean_equal <- function(outcome, group = NULL) {
  if (is.null(group)) {
    paste0("H0: the population mean or location of ", outcome, " equals the reference value.")
  } else {
    paste0("H0: the population mean or location of ", outcome, " is equal across levels of ", group, ".")
  }
}

h0_no_association <- function(x, y) {
  paste0("H0: ", x, " and ", y, " are independent.")
}

h0_no_correlation <- function(x, y) {
  paste0("H0: the correlation between ", x, " and ", y, " is 0.")
}

h0_proportion <- function(outcome, p) {
  paste0("H0: the true proportion of ", outcome, " successes equals ", p, ".")
}

formula_lhs_rhs <- function(expr, rhs_n = 1) {
  is_formula <- inherits(expr, "formula") || (is.call(expr) && identical(expr[[1]], as.name("~")))
  if (!is_formula || length(expr) != 3) {
    stop("Expected a two-sided formula such as outcome ~ group.", call. = FALSE)
  }
  lhs <- all.vars(expr[[2]])
  rhs <- all.vars(expr[[3]])
  if (length(lhs) != 1 || length(rhs) != rhs_n) {
    stop("Formula must contain one outcome and ", rhs_n, " right-hand side variable(s).", call. = FALSE)
  }
  list(lhs = lhs, rhs = rhs)
}

capture_pair <- function(outcome, group, rhs_label = "group") {
  expr <- substitute(outcome, parent.frame())
  if (inherits(expr, "formula")) {
    parsed <- formula_lhs_rhs(expr, rhs_n = 1)
    return(list(outcome = parsed$lhs, group = parsed$rhs))
  }
  if (missing(group)) {
    stop("Provide either `outcome` and `", rhs_label, "` columns or a formula like outcome ~ ", rhs_label, ".", call. = FALSE)
  }
  group_expr <- substitute(group, parent.frame())
  list(
    outcome = rlang::as_name(expr),
    group = rlang::as_name(group_expr)
  )
}

resolve_formula_pair <- function(first_expr, second_expr, third_expr, third_missing, rhs_label = "group") {
  if (inherits(first_expr, "formula") || (is.call(first_expr) && identical(first_expr[[1]], as.name("~")))) {
    parsed <- formula_lhs_rhs(first_expr, rhs_n = 1)
    return(list(outcome = parsed$lhs, group = parsed$rhs))
  }
  if (inherits(second_expr, "formula") || (is.call(second_expr) && identical(second_expr[[1]], as.name("~")))) {
    parsed <- formula_lhs_rhs(second_expr, rhs_n = 1)
    return(list(outcome = parsed$lhs, group = parsed$rhs))
  }
  if (third_missing) {
    stop("Provide a formula like outcome ~ ", rhs_label, " or provide separate outcome and ", rhs_label, " columns.", call. = FALSE)
  }
  list(
    outcome = rlang::as_name(second_expr),
    group = rlang::as_name(third_expr)
  )
}

resolve_data_first_or_formula <- function(first_value, second_value) {
  if (inherits(first_value, "data.frame")) first_value else second_value
}
