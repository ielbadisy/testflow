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
  meta <- workflow_meta(x)
  tf_with_cli_colors({
    tf_title(meta$title)
    tf_field(meta$outcome_label, x$outcome)
    tf_field(meta$group_label, x$group)
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
      tf_section("Effect size")
      if (is.na(x$effect_size$estimate[1])) {
        tf_bullet("Effect size not reported.")
      } else {
        tf_line(paste0(
          tf_label(x$effect_size$name[1]), " ",
          format_stat(x$effect_size$estimate[1]),
          if (!is.na(x$effect_size$magnitude[1])) paste0(", ", tf_value(x$effect_size$magnitude[1])) else ""
        ))
      }
      tf_blank()
    }

    if (!is.null(meta$table_data) && nrow(meta$table_data) > 0) {
      tf_section(meta$table_label)
      tf_table(meta$table_data)
      tf_blank()
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
  meta <- workflow_meta(x)
  tf_with_cli_colors({
    tf_title(paste0(meta$title, " (summary)"))
    print_summary_field("Workflow", x$workflow)
    print_summary_field(meta$outcome_label, x$outcome)
    print_summary_field(meta$group_label, x$group)
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

    if (!is.null(meta$table_data) && nrow(meta$table_data) > 0) {
      tf_blank()
      tf_section(meta$table_label)
      tf_table(meta$table_data)
    }

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

# Prints a small per-term/per-metric results table (regression
# coefficients, hazard ratios, the diagnostic/ROC/ICC/agreement tables,
# post-hoc comparisons, ...) to the console. Numeric columns are rounded
# with the same format_stat() 2-decimal convention used everywhere else in
# printed output; row names are dropped since none of these tables use them
# meaningfully. Uses base print.data.frame() rather than a new formatting
# dependency.
tf_table <- function(df) {
  df <- as.data.frame(df)
  p_cols <- grepl("^p(\\.value|\\.adj|value)?$", names(df))
  for (i in seq_along(df)) {
    if (is.integer(df[[i]])) {
      df[[i]] <- format(df[[i]], big.mark = ",", scientific = FALSE)
    } else if (is.numeric(df[[i]])) {
      df[[i]] <- if (p_cols[i]) format_p(df[[i]]) else format_stat(df[[i]])
    }
  }
  # Base data.frame printing wraps onto a second block of rows once the
  # formatted width exceeds getOption("width") (80 by default), which is
  # easy to trigger with a text column (e.g. long method/term labels)
  # alongside several numeric columns - and once it wraps, which
  # statistic/p-value belongs to which row is no longer visually obvious.
  # Widen the print width for just this call rather than truncating labels.
  col_widths <- mapply(function(col, nm) max(nchar(as.character(col)), nchar(nm), na.rm = TRUE), df, names(df))
  needed_width <- sum(col_widths) + 2 * ncol(df)
  old_width <- options(width = min(10000, max(getOption("width"), needed_width)))
  on.exit(options(old_width), add = TRUE)
  print(df, row.names = FALSE)
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
