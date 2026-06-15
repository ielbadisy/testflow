#' Print a testflow object
#' @param x A testflow object.
#' @param ... Unused.
#' @export
print.testflow <- function(x, ...) {
  cat("Statistical test workflow\n\n")
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
  h0 <- primary_h0(x)
  if (!is.na(h0)) cat(h0, "\n")
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
  ci <- primary_ci(x)
  ci_text <- if (any(is.na(ci))) "" else paste0(", 95% CI [", format_stat(ci[1]), ", ", format_stat(ci[2]), "]")
  paste0(
    "statistic = ", format_stat(stat),
    ifelse(is.na(df), "", paste0(", df = ", format_stat(df))),
    ", p = ", format_p(primary_p(x)),
    ci_text
  )
}

#' Summarize a testflow object
#' @param object A testflow object.
#' @param ... Unused.
#' @export
summary.testflow <- function(object, ...) {
  out <- list(
    question = paste(object$design, "workflow for", object$outcome %||% "outcome"),
    workflow = object$workflow,
    design = object$design,
    outcome = object$outcome,
    group = object$group,
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

#' @export
print.summary.testflow <- function(x, ...) {
  cat("testflow summary\n\n")
  print_summary_field("Workflow", x$workflow)
  print_summary_field("Design", x$design)
  print_summary_field("Outcome", x$outcome)
  print_summary_field("Group", x$group)
  print_summary_field("Recommended test", x$recommended$test %||% x$recommended)

  h0 <- x$primary_test$null_hypothesis %||% NA_character_
  print_summary_field("H0", h0)

  cat("\nResult\n")
  print_summary_field("Statistic", first_or_na(x$primary_test$statistic), formatter = format_stat)
  print_summary_field("df", first_or_na(x$primary_test$parameter), formatter = format_stat)
  print_summary_field("p", first_or_na(x$primary_test$p.value), formatter = format_p)

  if (!is.null(x$primary_test$conf.low) && !is.null(x$primary_test$conf.high)) {
    ci <- paste0("[", format_stat(x$primary_test$conf.low[[1]]), ", ", format_stat(x$primary_test$conf.high[[1]]), "]")
    print_summary_field("95% CI", ci)
  }

  if (!is.null(x$effect_size) && nrow(x$effect_size) > 0) {
    effect <- paste0(
      x$effect_size$name[1], " = ", format_stat(x$effect_size$estimate[1]),
      if (!is.na(x$effect_size$magnitude[1])) paste0(", ", x$effect_size$magnitude[1]) else ""
    )
    print_summary_field("Effect size", effect)
  }

  print_summary_field("Decision", x$decision)

  if (!is.null(x$report) && !is.na(x$report)) {
    cat("\nReport\n")
    cat(x$report, "\n")
  }

  invisible(x)
}

print_summary_field <- function(label, value, formatter = identity) {
  value <- first_or_na(value)
  if (length(value) == 0 || is.null(value) || is.na(value)) return(invisible(NULL))
  cat(sprintf("%-18s %s\n", paste0(label, ":"), formatter(value)))
  invisible(NULL)
}

first_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA)
  unname(x[[1]])
}
