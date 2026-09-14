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
  tibble::as_tibble(basetable::rbindfill(basetable::map(.x, .f, ...), fill = TRUE))
}

.imap_dfr <- function(.x, .f) {
  idx <- if (!is.null(names(.x))) names(.x) else seq_along(.x)
  tibble::as_tibble(basetable::rbindfill(basetable::traverse(list(.x, idx), .f), fill = TRUE))
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
