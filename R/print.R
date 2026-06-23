#' Print a testflow object
#'
#' @description
#' Print a testflow object.
#'
#' @details
#' Console colors are enabled by default in interactive sessions. Use
#' `options(testflow.cli_colors = FALSE)` to disable colors, or
#' `options(testflow.cli_colors = TRUE)` to force colors in non-interactive
#' output.
#'
#' @param x A testflow object.
#' @param ... Unused.
#' @return The input `testflow` object, invisibly. Called for its side effect of
#' printing a formatted workflow summary to the console.
#' @export
print.testflow <- function(x, ...) {
  tf_with_cli_colors({
    tf_title("Statistical test workflow")
    tf_field("Outcome", x$outcome)
    tf_field("Group", x$group)
    tf_field("Design", x$design)
    tf_blank()

    if (!is.null(x$assumptions)) {
      tf_section("Assumptions")
      print_assumption_line(x$assumptions)
      tf_blank()
    }

    tf_section("Recommended test")
    tf_line(tf_value(x$recommended$test %||% x$recommended))
    tf_blank()

    tf_section("Result")
    h0 <- primary_h0(x)
    if (!is.na(h0)) tf_line(h0)
    tf_line(format_primary_result(x))
    tf_blank()

    if (!is.null(x$effect_size) && nrow(x$effect_size) > 0) {
      if (is.na(x$effect_size$estimate[1]) || is.na(x$effect_size$magnitude[1])) {
        tf_section("Effect size")
        tf_bullet("Effect size not reported.")
        tf_blank()
      } else {
      tf_section("Effect size")
      tf_line(paste0(
        tf_label(x$effect_size$name[1]), " ",
        format_stat(x$effect_size$estimate[1]), ", ",
        tf_value(x$effect_size$magnitude[1])
      ))
      tf_blank()
      }
    }

    tf_section("Report")
    tf_line(report(x))
    invisible(x)
  })
}

print_assumption_line <- function(assumptions) {
  assumptions <- format_assumptions(assumptions)
  if (!nrow(assumptions)) {
    tf_bullet("No assumptions reported.")
    return(invisible(NULL))
  }
  if (!"name" %in% names(assumptions)) {
    assumptions$name <- paste0("Assumption ", seq_len(nrow(assumptions)))
  } else {
    assumptions$name <- ifelse(is.na(assumptions$name) | !nzchar(as.character(assumptions$name)), paste0("Assumption ", seq_len(nrow(assumptions))), assumptions$name)
  }
  assumptions <- dplyr::select(assumptions, dplyr::any_of(c("name", "status", "message", "method", "statistic", "p_value", "details")))
  purrr::pwalk(assumptions, function(name, status, message, method, statistic, p_value, details) {
    line <- paste0(tf_label(name), " ", tf_value(status), ": ", message)
    extras <- c()
    if (!is.na(method) && nzchar(method)) extras <- c(extras, paste0("method=", method))
    if (!is.na(statistic)) extras <- c(extras, paste0("statistic=", format_stat(statistic)))
    if (!is.na(p_value)) extras <- c(extras, paste0("p=", format_p(p_value)))
    if (!is.na(details) && nzchar(details)) extras <- c(extras, details)
    if (length(extras)) line <- paste0(line, " (", paste(extras, collapse = "; "), ")")
    tf_bullet(line)
  })
  invisible(NULL)
}

format_primary_result <- function(x) {
  if (identical(x$workflow, "correlation_matrix")) {
    return(paste0("pairwise correlations reported; smallest pairwise p = ", format_p(primary_p(x))))
  }
  if (identical(x$workflow, "outliers")) {
    return(paste0("flagged rows = ", format_count(primary_statistic(x))))
  }

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
#'
#' @description
#' Summarize a testflow object.
#'
#' @details
#' Console colors follow the same `testflow.cli_colors` option used by
#' [print.testflow()].
#'
#' @param object A testflow object.
#' @param ... Unused.
#' @return A `summary.testflow` list containing the workflow metadata,
#' descriptives, assumptions, recommended test, primary and alternative test
#' results, post-hoc results when available, effect size, decision, and report
#' text.
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
  tf_with_cli_colors({
    tf_title("testflow summary")
    print_summary_field("Workflow", x$workflow)
    print_summary_field("Design", x$design)
    print_summary_field("Outcome", x$outcome)
    print_summary_field("Group", x$group)
    print_summary_field("Recommended test", x$recommended$test %||% x$recommended)

    h0 <- x$primary_test$null_hypothesis %||% NA_character_
    print_summary_field("H0", h0)

    tf_blank()
    tf_section("Result")
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
      tf_blank()
      tf_section("Report")
      tf_line(x$report)
    }

    invisible(x)
  })
}

print_summary_field <- function(label, value, formatter = identity) {
  value <- first_or_na(value)
  if (length(value) == 0 || is.null(value) || is.na(value)) return(invisible(NULL))
  tf_line(paste0(tf_label(label), " ", tf_value(formatter(value))))
  invisible(NULL)
}

first_or_na <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA)
  unname(x[[1]])
}

format_count <- function(x) {
  if (is.na(x)) return("NA")
  format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
}

tf_with_cli_colors <- function(expr) {
  color_setting <- getOption("testflow.cli_colors", NULL)
  old_options <- NULL

  if (isTRUE(color_setting) || (is.null(color_setting) && interactive())) {
    old_options <- options(cli.num_colors = 256)
  } else if (isFALSE(color_setting)) {
    old_options <- options(cli.num_colors = 1)
  }

  if (!is.null(old_options)) {
    on.exit(options(old_options), add = TRUE)
  }

  eval.parent(substitute(expr))
}

tf_title <- function(text) {
  tf_line(cli::style_bold(cli::col_blue(text)))
  tf_blank()
}

tf_section <- function(text) {
  tf_line(cli::style_bold(cli::col_magenta(text)))
}

tf_line <- function(...) {
  cat(..., "\n", sep = "")
}

tf_blank <- function() {
  cat("\n")
}

tf_bullet <- function(text) {
  tf_line("* ", text)
}

tf_field <- function(label, value) {
  if (length(value) == 0 || is.null(value) || is.na(value)) return(invisible(NULL))
  tf_line(paste0(tf_label(label), " ", tf_value(value)))
  invisible(NULL)
}

tf_label <- function(label) {
  cli::col_cyan(paste0(label, ":"))
}

tf_value <- function(value) {
  cli::col_green(as.character(value))
}
