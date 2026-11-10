#' Convert a testflow result to a one-row table
#'
#' `as.data.frame()` returns the one-row summary of a testflow or sample-size
#' result as a plain data frame. `as_tibble()` is kept for existing code and
#' returns the same row as a basetable; testflow no longer depends on
#' 'tibble'.
#'
#' @param x A testflow object or a sample-size object.
#' @param row.names,optional Unused; present for compatibility with
#'   [base::as.data.frame()].
#' @param ... Unused.
#' @return A one-row table with the workflow name, design, variables,
#' recommended test, null hypothesis, statistic, degrees of freedom when
#' available, p-value, confidence interval when available, effect-size fields,
#' and decision text. For a sample-size object: endpoint, design, objective,
#' method, sample sizes, report and any diagnostics.
#' @name as_tibble
#' @export
as_tibble <- function(x, ...) {
  UseMethod("as_tibble")
}

#' @rdname as_tibble
#' @export
as.data.frame.testflow <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- as_tibble.testflow(x)
  class(out) <- "data.frame"
  out
}

#' @rdname as_tibble
#' @export
as.data.frame.sample_size <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- as_tibble.sample_size(x)
  class(out) <- "data.frame"
  out
}

#' @rdname as_tibble
#' @export
as_tibble.testflow <- function(x, ...) {
  effect <- x$effect_size
  ci <- primary_ci(x)
  .bt(
    workflow = x$workflow,
    design = x$design,
    outcome = x$outcome %||% NA_character_,
    group = x$group %||% NA_character_,
    recommended_test = x$recommended$test %||% as.character(x$recommended),
    null_hypothesis = primary_h0(x),
    statistic = primary_statistic(x),
    df = primary_df(x),
    p = primary_p(x),
    conf.low = ci[1],
    conf.high = ci[2],
    effect_size = if (!is.null(effect) && nrow(effect) > 0) effect$estimate[1] else NA_real_,
    effect_size_name = if (!is.null(effect) && nrow(effect) > 0) effect$name[1] else NA_character_,
    effect_size_magnitude = if (!is.null(effect) && nrow(effect) > 0) effect$magnitude[1] else NA_character_,
    decision = significance_decision(primary_p(x))
  )
}
