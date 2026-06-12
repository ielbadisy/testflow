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
