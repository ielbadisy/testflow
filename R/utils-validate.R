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
  if (!"statistic" %in% names(out)) out$statistic <- NA_real_
  if (!"parameter" %in% names(out)) out$parameter <- NA_real_
  if (!is.null(method)) {
    out$method <- method
  }
  out
}
