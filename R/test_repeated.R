#' Run a repeated-measures workflow from wide data
#' @param data A data frame.
#' @param measures Repeated numeric columns selected with tidyselect syntax.
#' @param id Optional subject identifier.
#' @param between Optional between-subject factor (constant within each
#'   subject). When supplied, a mixed (split-plot) ANOVA,
#'   `outcome ~ within * between + Error(id / within)`, is added as
#'   `$alternative_tests$mixed_anova` (one row per term: between, within, and
#'   their interaction) and summarized in the report. The primary result stays
#'   the within-subject test (repeated-measures ANOVA or Friedman).
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_repeated`. The object is a
#' list containing long-form repeated-measures data, numeric descriptives,
#' assumption checks, recommended repeated-measures ANOVA or Friedman test,
#' primary test result with null hypothesis, alternative repeated-measures
#' results, post-hoc paired comparisons, effect size, optional `ggplot`,
#' original call, and report text.
#' @references
#' Fisher, R. A. (1925). \emph{Statistical Methods for Research Workers}.
#' Oliver and Boyd.
#'
#' Friedman, M. (1937). The use of ranks to avoid the assumption of normality
#' implicit in the analysis of variance. \emph{Journal of the American
#' Statistical Association}, 32(200), 675-701.
#'
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods.
#' \emph{Biometrics Bulletin}, 1(6), 80-83.
#'
#' Girden, E. R. (1992). \emph{ANOVA: Repeated Measures}. Sage.
#' @export
testrepeated <- function(data, measures, id = NULL, between = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  measure_nms <- tidyselect_names(data, {{ measures }})
  warn_if_two_measures_repeated(measure_nms)
  id_nm <- if (missing(id) || is.null(substitute(id))) ".testflow_id" else rlang::as_name(rlang::ensym(id))
  df <- data
  if (id_nm == ".testflow_id") df[[id_nm]] <- seq_len(nrow(df))
  long <- .pivot_longer_simple(df, measure_nms, "time", "value")
  between_nm <- if (is.null(substitute(between))) NULL else rlang::as_name(rlang::ensym(between))
  repeated_core(long, "value", "time", id_nm, outcome_label = paste(measure_nms, collapse = ", "), display_by_within = TRUE, alpha = alpha, plot = plot, na.rm = na.rm, call = match.call(), between_nm = between_nm)
}

#' Run a repeated-measures workflow from long data
#' @param data A data frame.
#' @param outcome Numeric outcome column.
#' @param within Within-subject time/condition column.
#' @param id Subject identifier column.
#' @param between Optional between-subject factor (constant within each
#'   subject). When supplied, a mixed (split-plot) ANOVA,
#'   `outcome ~ within * between + Error(id / within)`, is added as
#'   `$alternative_tests$mixed_anova` (one row per term: between, within, and
#'   their interaction) and summarized in the report. The primary result stays
#'   the within-subject test (repeated-measures ANOVA or Friedman).
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_repeated`. The object is a
#' list containing the cleaned long-format data, numeric descriptives,
#' assumption checks, recommended repeated-measures ANOVA or Friedman test,
#' primary test result with null hypothesis, alternative repeated-measures
#' results, post-hoc paired comparisons, effect size, optional `ggplot`,
#' original call, and report text.
#' @export
testrepeatedlong <- function(data, outcome, within, id, between = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  within_nm <- rlang::as_name(rlang::ensym(within))
  id_nm <- rlang::as_name(rlang::ensym(id))
  between_nm <- if (is.null(substitute(between))) NULL else rlang::as_name(rlang::ensym(between))
  repeated_core(data, outcome_nm, within_nm, id_nm, outcome_label = outcome_nm, display_by_within = FALSE, alpha = alpha, plot = plot, na.rm = na.rm, call = match.call(), between_nm = between_nm)
}

repeated_core <- function(data, outcome_nm, within_nm, id_nm, outcome_label = NULL, display_by_within = FALSE, alpha = 0.05, plot = TRUE, na.rm = TRUE, call = NULL, between_nm = NULL) {
  if (is.null(outcome_label)) outcome_label <- outcome_nm
  df <- drop_missing(data, c(outcome_nm, within_nm, id_nm, between_nm), na.rm = na.rm)
  # aov()'s Error(id/within) requires id to be a factor to build the correct
  # id:within error stratum; a numeric/integer id (e.g. the auto-generated
  # subject index) silently produces the wrong stratum and wrong error df.
  df[[id_nm]] <- as.factor(df[[id_nm]])
  normality <- check_normality(df, outcome_nm, within_nm, alpha)
  if (isTRUE(display_by_within)) {
    normality$name <- sub(paste0("^Normality: ", outcome_nm, " \\((.*)\\)$"), "Normality: \\1", normality$name)
    normality$variable <- normality$group
  }
  wide <- stats::reshape(
    as.data.frame(df[, c(id_nm, within_nm, outcome_nm), drop = FALSE]),
    idvar = id_nm, timevar = within_nm, direction = "wide"
  )
  # stats::reshape() names the new value columns "<outcome_nm>.<level>"; strip
  # the fixed prefix to match tidyr::pivot_wider()'s bare-level column names.
  value_cols <- setdiff(names(wide), id_nm)
  names(wide)[match(value_cols, names(wide))] <- substring(value_cols, nchar(outcome_nm) + 2)
  rownames(wide) <- NULL
  wide <- .as_bt(wide)
  measure_nms <- setdiff(names(wide), id_nm)
  complete <- stats::na.omit(wide[, measure_nms, drop = FALSE])
  sphericity <- check_sphericity(as.matrix(complete), alpha)
  friedman <- stats::friedman.test(as.matrix(complete))
  anova <- repeated_anova_test(df, outcome_nm, within_nm, id_nm)
  recommendation <- if (all(normality$status == "acceptable")) "Repeated-measures ANOVA" else "Friedman test"
  h0 <- h0_mean_equal(outcome_label, within_nm)
  primary <- if (recommendation == "Friedman test") {
    add_null_hypothesis(safe_tidy_htest(friedman, "Friedman test"), h0)
  } else {
    add_null_hypothesis(anova, h0)
  }
  effect <- if (recommendation == "Friedman test") {
    friedman_effect_size(friedman, n = nrow(complete), k = length(measure_nms))
  } else {
    repeated_eta_squared(df, outcome_nm, within_nm, id_nm)
  }
  posthoc_method <- if (recommendation == "Friedman test") "wilcox" else "t"
  posthoc_tests <- paired_posthoc_numeric(complete, method = posthoc_method, p.adjust.method = "BH")
  plt <- if (plot) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[within_nm]], y = .data[[outcome_nm]], group = .data[[id_nm]])) +
      ggplot2::geom_line(alpha = 0.25) +
      ggplot2::geom_point(alpha = 0.45) +
      ggplot2::stat_summary(ggplot2::aes(group = 1), fun = mean, geom = "line", linewidth = 1.1, color = "#F58518") +
      ggplot2::labs(title = "Repeated-measures workflow", subtitle = paste0(recommendation, ", p = ", format_p(primary$p.value[1])), x = within_nm, y = outcome_nm) +
      ggplot2::theme_minimal()
  } else NULL
  alt <- list(anova = add_null_hypothesis(anova, h0), friedman = add_null_hypothesis(safe_tidy_htest(friedman, "Friedman test"), h0))
  if (!is.null(between_nm)) alt$mixed_anova <- mixed_anova_test(df, outcome_nm, within_nm, id_nm, between_nm)
  out <- new_testflow("repeated", "repeated numeric measurements", outcome_label, within_nm, id_nm, df, descriptives_numeric(df, outcome_nm, within_nm), assumption_checks(check_independence_note("Repeated measurements from the same subjects are assumed by design."), normality, sphericity), list(test = recommendation), primary, alt, posthoc = posthoc_tests, effect_size = effect, plot = plt, call = call, subclass = "repeated")
  out$interpretation <- make_report(out, alpha)
  if (!is.null(alt$mixed_anova)) {
    m <- alt$mixed_anova
    out$interpretation <- paste0(
      out$interpretation, " Mixed (split-plot) ANOVA with between-subject factor ", between_nm, ": ",
      paste0(m$term, ": F = ", format_stat(m$statistic), ", df = ", format_stat(m$parameter), ", ", format_stat(m$den.df), ", p = ", format_p(m$p.value), collapse = "; "),
      "."
    )
  }
  out
}

# Mixed (split-plot) ANOVA: outcome ~ within * between + Error(id / within).
# The between-subject term is tested in the subject stratum; the within term
# and the interaction in the subject-by-within stratum.
mixed_anova_test <- function(long, outcome, within, id, between) {
  long[[between]] <- as.factor(long[[between]])
  per_id <- tapply(as.character(long[[between]]), long[[id]], function(v) length(unique(v)))
  if (any(per_id > 1, na.rm = TRUE)) {
    stop("`between` must be constant within each subject.", call. = FALSE)
  }
  formula <- stats::as.formula(paste(
    tf_backtick(outcome), "~", tf_backtick(within), "*", tf_backtick(between),
    "+ Error(", tf_backtick(id), "/", tf_backtick(within), ")"
  ))
  fit <- stats::aov(formula, data = long)
  strata <- summary(fit)
  rows <- lapply(strata, function(st) {
    tb <- st[[1]]
    terms <- trimws(rownames(tb))
    is_res <- terms == "Residuals"
    if (!any(is_res) || !any(!is_res)) return(NULL)
    keep <- !is_res & !is.na(tb[["F value"]])
    .bt(
      method = "Mixed (split-plot) ANOVA",
      term = terms[keep],
      statistic = tb[["F value"]][keep],
      parameter = tb[["Df"]][keep],
      den.df = rep(tb[["Df"]][is_res][1], sum(keep)),
      p.value = tb[["Pr(>F)"]][keep]
    )
  })
  out <- .rbind_rows(rows)
  out$null_hypothesis <- ifelse(
    grepl(":", out$term, fixed = TRUE),
    paste0("H0: the change in ", outcome, " across ", within, " is the same at every level of ", between, " (no interaction)."),
    paste0("H0: the mean of ", outcome, " is equal across levels of ", out$term, ".")
  )
  out
}
