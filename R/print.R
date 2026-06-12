#' Print a testflow object
#' @param x A testflow object.
#' @param ... Unused.
#' @export
print.testflow <- function(x, ...) {
  cat("Clinical statistical test workflow\n\n")
  if (!is.null(x$outcome)) cat("Outcome:", x$outcome, "\n")
  if (!is.null(x$group)) cat("Group:", x$group, "\n")
  cat("Design:", x$design, "\n\n")

  if (!is.null(x$assumptions)) {
    cat("Assumptions:\n")
    print_assumption_line(x$assumptions)
    cat("\n")
  }

  cat("Recommended test:\n")
  cat(x$recommended$test %||% x$recommended, "\n\n")

  cat("Result:\n")
  cat(format_primary_result(x), "\n\n")

  if (!is.null(x$effect_size) && nrow(x$effect_size) > 0) {
    cat("Effect size:\n")
    cat(x$effect_size$name[1], " = ", format_stat(x$effect_size$estimate[1]), ", ", x$effect_size$magnitude[1], "\n\n", sep = "")
  }

  cat("Report:\n")
  cat(report(x), "\n")
  invisible(x)
}

print_assumption_line <- function(assumptions) {
  if (is.list(assumptions) && !inherits(assumptions, "data.frame")) {
    for (nm in names(assumptions)) {
      value <- assumptions[[nm]]
      if (inherits(value, "data.frame") && "status" %in% names(value)) {
        cat("- ", nm, ": ", paste(unique(value$status), collapse = ", "), "\n", sep = "")
      }
    }
  } else if (inherits(assumptions, "data.frame") && "status" %in% names(assumptions)) {
    cat("- Status: ", paste(unique(assumptions$status), collapse = ", "), "\n", sep = "")
  }
}

format_primary_result <- function(x) {
  test <- x$primary_test
  if (is.null(test)) return("No primary test available.")
  stat <- primary_statistic(x)
  df <- primary_df(x)
  paste0(
    "statistic = ", format_stat(stat),
    ifelse(is.na(df), "", paste0(", df = ", format_stat(df))),
    ", p = ", format_p(primary_p(x))
  )
}

#' Summarize a testflow object
#' @param object A testflow object.
#' @param ... Unused.
#' @export
summary.testflow <- function(object, ...) {
  out <- list(
    question = paste(object$design, "workflow for", object$outcome %||% "outcome"),
    design = object$design,
    descriptives = object$descriptives,
    assumptions = object$assumptions,
    recommended = object$recommended,
    primary_test = object$primary_test,
    alternative_tests = object$alternative_tests,
    posthoc = object$posthoc,
    effect_size = object$effect_size,
    decision = significance_decision(primary_p(object)),
    report = report(object)
  )
  class(out) <- "summary.testflow"
  out
}
