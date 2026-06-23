#' Convert a testflow object to a one-row tibble
#' @param x A testflow object.
#' @param ... Unused.
#' @return A one-row tibble with the workflow name, design, variables,
#' recommended test, null hypothesis, statistic, degrees of freedom when
#' available, p-value, confidence interval when available, effect-size fields,
#' and decision text.
#' @export
as_tibble <- function(x, ...) {
  UseMethod("as_tibble")
}

#' @export
as_tibble.testflow <- function(x, ...) {
  effect <- x$effect_size
  ci <- primary_ci(x)
  tibble::tibble(
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
