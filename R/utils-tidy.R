new_testflow <- function(
  workflow,
  design,
  outcome = NULL,
  group = NULL,
  id = NULL,
  data = NULL,
  descriptives = NULL,
  assumptions = NULL,
  recommended = NULL,
  primary_test = NULL,
  alternative_tests = NULL,
  posthoc = NULL,
  effect_size = NULL,
  interpretation = NULL,
  plot = NULL,
  call = NULL,
  subclass = NULL
) {
  x <- list(
    workflow = workflow,
    design = design,
    outcome = outcome,
    group = group,
    id = id,
    data = data,
    descriptives = descriptives,
    assumptions = assumptions,
    recommended = recommended,
    primary_test = add_effect_size_result(primary_test, effect_size),
    alternative_tests = alternative_tests,
    posthoc = posthoc,
    effect_size = effect_size,
    interpretation = interpretation,
    plot = plot,
    call = call
  )

  class(x) <- c(paste0("testflow_", subclass %||% workflow), "testflow")
  x
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# basetable/base-R replacements for purrr::map_dfr()/purrr::imap_dfr(): map
# (or index-map) .f over .x and row-bind the results, matching purrr's
# silent drop of NULL results (basetable::rbindfill() already drops NULLs).
.map_dfr <- function(.x, .f, ...) {
  .rbind_rows(basetable::map(.x, .f, ...))
}

.imap_dfr <- function(.x, .f) {
  idx <- if (!is.null(names(.x))) names(.x) else seq_along(.x)
  .rbind_rows(basetable::traverse(list(.x, idx), .f))
}

# Row-bind data frames like dplyr::bind_rows(). basetable::rbindfill() stops on
# any column type mismatch; before binding, widen only the safe cases:
# integer + double -> double (e.g. median() is integer for an odd-sized group
# and double for an even-sized one), and all-NA logical columns -> the type of
# the other inputs. Other mismatches (e.g. numeric vs character) still error.
.rbind_rows <- function(...) {
  dfs <- list(...)
  if (length(dfs) == 1L && is.list(dfs[[1]]) && !is.data.frame(dfs[[1]])) dfs <- dfs[[1]]
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0L) return(tibble::tibble())
  cols <- unique(unlist(lapply(dfs, names)))
  for (col in cols) {
    has <- vapply(dfs, function(d) col %in% names(d), logical(1))
    vals <- lapply(dfs[has], `[[`, col)
    all_na_lgl <- vapply(vals, function(v) is.logical(v) && all(is.na(v)), logical(1))
    types <- unique(vapply(vals[!all_na_lgl], typeof, character(1)))
    target <- if (length(types) == 0L) {
      NULL
    } else if (all(types %in% c("integer", "double")) && length(types) == 2L) {
      "double"
    } else if (length(types) == 1L && any(all_na_lgl)) {
      types
    } else {
      NULL
    }
    if (is.null(target)) next
    idx <- which(has)
    for (i in idx) {
      v <- dfs[[i]][[col]]
      if (is.factor(v) || inherits(v, c("Date", "POSIXt"))) next
      if (target == "double" && is.integer(v)) {
        dfs[[i]][[col]] <- as.double(v)
      } else if (is.logical(v) && all(is.na(v)) && target != "logical") {
        dfs[[i]][[col]] <- as.vector(v, mode = target)
      }
    }
  }
  tibble::as_tibble(basetable::rbindfill(dfs, fill = TRUE))
}

# tidyr::pivot_longer() replacement for the simple "melt these columns, keep
# everything else as id columns" case used throughout this package (no
# multi-value-column / names_pattern features).
.pivot_longer_simple <- function(df, cols, names_to, values_to) {
  id_cols <- setdiff(names(df), cols)
  n <- nrow(df)
  k <- length(cols)
  row_idx <- rep(seq_len(n), each = k)
  id_part <- df[row_idx, id_cols, drop = FALSE]
  name_part <- rep(cols, times = n)
  value_cols <- df[, cols, drop = FALSE]
  value_part <- do.call(c, lapply(seq_len(n), function(i) unlist(value_cols[i, ], use.names = FALSE)))
  out <- cbind(id_part, stats::setNames(
    data.frame(name_part, value_part, stringsAsFactors = FALSE),
    c(names_to, values_to)
  ))
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

primary_p <- function(x) {
  if (is.null(x$primary_test$p.value)) NA_real_ else x$primary_test$p.value
}

primary_statistic <- function(x) {
  if (is.null(x$primary_test$statistic)) NA_real_ else unname(x$primary_test$statistic[[1]])
}

primary_df <- function(x) {
  if (!is.null(x$primary_test$parameter)) unname(x$primary_test$parameter[[1]]) else NA_real_
}

primary_ci <- function(x) {
  if (is.null(x$primary_test$conf.low) || is.null(x$primary_test$conf.high)) {
    return(c(NA_real_, NA_real_))
  }
  c(x$primary_test$conf.low[[1]], x$primary_test$conf.high[[1]])
}

primary_h0 <- function(x) {
  if (is.null(x$primary_test$null_hypothesis)) NA_character_ else x$primary_test$null_hypothesis[[1]]
}

# Fisher's exact test that does not crash on large tables: when the exact
# network algorithm runs out of workspace (FEXACT error), fall back to a
# Monte Carlo p-value (B = 10000) and say so in the method label.
.fisher_safe <- function(tab) {
  res <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
  if (!is.null(res)) return(res)
  res <- stats::fisher.test(tab, simulate.p.value = TRUE, B = 10000)
  res$method <- "Fisher exact test (Monte Carlo p-value, B = 10000)"
  res
}

# Shapiro-Wilk on model residuals with the same guards as check_normality():
# "not enough data" below 3 observations, "not checked" above 5000
# (stats::shapiro.test() errors outside 3..5000).
.residual_normality <- function(resid, alpha, ok_text, bad_text) {
  resid <- resid[!is.na(resid)]
  n <- length(resid)
  if (n < 3 || n > 5000) {
    if (n > 5000) warning("Shapiro-Wilk is not recommended for n > 5000.", call. = FALSE)
    return(assumption_check(
      "Normality of residuals",
      if (n < 3) "not enough data" else "not checked",
      if (n < 3) "Too few observations for a useful normality check." else "Shapiro-Wilk is not recommended for n > 5000; inspect the residual plot instead.",
      method = "Shapiro-Wilk", statistic = NA_real_, p_value = NA_real_
    ))
  }
  test <- stats::shapiro.test(resid)
  assumption_check(
    "Normality of residuals",
    ifelse(test$p.value >= alpha, "acceptable", "not acceptable"),
    ifelse(test$p.value >= alpha, ok_text, bad_text),
    method = "Shapiro-Wilk",
    statistic = unname(test$statistic),
    p_value = test$p.value
  )
}
