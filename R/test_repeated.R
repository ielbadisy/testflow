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
#' @param times Optional labels for the measures on the plot x axis, one per
#'   measure in order (for example `0:4` months); the axis is then titled
#'   "Time". Only the plot changes.
#' @param xlab,ylab,legendtitle Optional plot labels for the x axis, the y
#'   axis and the legend of the `between` groups. Left out, the defaults are
#'   kept; a character value replaces them; `NULL` removes them.
#' @param posthoc Post-hoc method. `"auto"` (default) keeps paired t tests
#'   after the repeated-measures ANOVA and paired Wilcoxon signed-rank tests
#'   after the Friedman test. `"t"` or `"wilcoxon"` fixes them; `"nemenyi"`
#'   and `"conover"` give the all-pairs tests on within-subject ranks that
#'   follow a Friedman test (see Details); `"none"` skips post hoc tests. All
#'   use the subjects with every measure.
#' @param p.adjust.method Multiplicity adjustment passed to [stats::p.adjust()]
#'   for the `"t"`, `"wilcoxon"` and `"conover"` post hoc tests (Nemenyi
#'   carries its own studentized-range adjustment).
#' @return A `testflow` object with class `testflow_repeated`. The object is a
#' list containing long-form repeated-measures data, numeric descriptives,
#' assumption checks, recommended repeated-measures ANOVA or Friedman test,
#' primary test result with null hypothesis, alternative repeated-measures
#' results, post-hoc paired comparisons, effect size, optional `ggplot`,
#' original call, and report text.
#'
#' The repeated-measures ANOVA rows (`$alternative_tests$anova`, and
#' `$primary_test` when it is the ANOVA) carry `gg.eps` and `hf.eps`
#' (Greenhouse-Geisser and Huynh-Feldt epsilon, from the subjects with all
#' measures) and `p.gg`, `p.hf`, the p-values with both df multiplied by
#' epsilon (Huynh-Feldt epsilon capped at 1), as `anova.mlm(test =
#' "Spherical")`. The mixed ANOVA rows carry the same columns for the within
#' term and the interaction (epsilon from the pooled within-group covariance;
#' `NA` for the between term).
#' @details
#' Rank post hoc tests. Each subject's \eqn{k} measures are ranked (ties get
#' mid-ranks); \eqn{R_j} is the rank sum and \eqn{\bar R_j = R_j / n} the mean
#' rank of condition \eqn{j} over the \eqn{n} complete subjects. `estimate` is
#' \eqn{\bar R_j - \bar R_i} (group2 minus group1).
#'
#' Nemenyi (1963): \deqn{q_{ij} = \sqrt{2}\,\frac{\bar R_j - \bar R_i}{\sqrt{k(k+1)/(6n)}},}
#' with p from the studentized range with \eqn{k} means and infinite df
#' (single step, `p.adj = p`). The variance is the untied one, so with many
#' within-subject ties the test is conservative.
#'
#' Conover (1999): with \eqn{A = \sum r^2} over all ranks,
#' \deqn{t_{ij} = \frac{R_j - R_i}{\sqrt{2 (n A - \sum_j R_j^2) / ((n-1)(k-1))}},}
#' referred to a t distribution with \eqn{(n-1)(k-1)} df, then adjusted by
#' `p.adjust.method`. \eqn{A} is computed from the observed ranks, which
#' corrects for ties.
#'
#' Both match `PMCMRplus::frdAllPairsNemenyiTest()` and
#' `PMCMRplus::frdAllPairsConoverTest()`.
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
#'
#' Greenhouse, S. W., and Geisser, S. (1959). On methods in the analysis of
#' profile data. \emph{Psychometrika}, 24(2), 95-112.
#'
#' Huynh, H., and Feldt, L. S. (1976). Estimation of the Box correction for
#' degrees of freedom from sample data in randomized block and split-plot
#' designs. \emph{Journal of Educational Statistics}, 1(1), 69-82.
#'
#' Nemenyi, P. (1963). \emph{Distribution-free Multiple Comparisons}. PhD
#' thesis, Princeton University.
#'
#' Conover, W. J. (1999). \emph{Practical Nonparametric Statistics} (3rd ed.).
#' Wiley.
#' @export
testrepeated <- function(data, measures, id = NULL, between = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE,
                         times = NULL, xlab, ylab, legendtitle,
                         posthoc = c("auto", "t", "wilcoxon", "nemenyi", "conover", "none"),
                         p.adjust.method = "BH") {
  measure_nms <- tidyselect_names(data, {{ measures }})
  warn_if_two_measures_repeated(measure_nms)
  id_nm <- if (missing(id) || is.null(substitute(id))) ".testflow_id" else rlang::as_name(rlang::ensym(id))
  df <- data
  if (id_nm == ".testflow_id") df[[id_nm]] <- seq_len(nrow(df))
  long <- .pivot_longer_simple(df, measure_nms, "time", "value")
  between_nm <- if (is.null(substitute(between))) NULL else rlang::as_name(rlang::ensym(between))
  if (!is.null(times) && length(times) != length(measure_nms)) stop("`times` needs one label per measure.", call. = FALSE)
  out <- repeated_core(long, "value", "time", id_nm, outcome_label = paste(measure_nms, collapse = ", "), display_by_within = TRUE, alpha = alpha, plot = plot, na.rm = na.rm, call = match.call(), between_nm = between_nm,
                       posthoc = match.arg(posthoc), p.adjust.method = p.adjust.method)
  if (!is.null(times) && !is.null(out$plot))
    out$plot <- out$plot + ggplot2::scale_x_discrete(labels = stats::setNames(as.character(times), measure_nms)) +
      ggplot2::labs(x = "Time")
  out$plot <- relabel(out$plot, if (!missing(xlab)) list(x = xlab), if (!missing(ylab)) list(y = ylab),
                      if (!missing(legendtitle)) list(colour = legendtitle, shape = legendtitle))
  out
}

# Apply only the label overrides given (an explicit NULL removes the label).
relabel <- function(p, ...) {
  lab <- c(...)
  if (is.null(p) || !length(lab)) p else p + do.call(ggplot2::labs, lab)
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
#' @param xlab,ylab,legendtitle Optional plot labels for the x axis, the y
#'   axis and the legend of the `between` groups. Left out, the defaults are
#'   kept; a character value replaces them; `NULL` removes them.
#' @param posthoc Post-hoc method. `"auto"` (default) keeps paired t tests
#'   after the repeated-measures ANOVA and paired Wilcoxon signed-rank tests
#'   after the Friedman test. `"t"` or `"wilcoxon"` fixes them; `"nemenyi"`
#'   and `"conover"` give the all-pairs tests on within-subject ranks that
#'   follow a Friedman test (see Details); `"none"` skips post hoc tests. All
#'   use the subjects with every measure.
#' @param p.adjust.method Multiplicity adjustment passed to [stats::p.adjust()]
#'   for the `"t"`, `"wilcoxon"` and `"conover"` post hoc tests (Nemenyi
#'   carries its own studentized-range adjustment).
#' @return A `testflow` object with class `testflow_repeated`. The object is a
#' list containing the cleaned long-format data, numeric descriptives,
#' assumption checks, recommended repeated-measures ANOVA or Friedman test,
#' primary test result with null hypothesis, alternative repeated-measures
#' results, post-hoc paired comparisons, effect size, optional `ggplot`,
#' original call, and report text.
#'
#' The repeated-measures ANOVA rows (`$alternative_tests$anova`, and
#' `$primary_test` when it is the ANOVA) carry `gg.eps` and `hf.eps`
#' (Greenhouse-Geisser and Huynh-Feldt epsilon, from the subjects with all
#' measures) and `p.gg`, `p.hf`, the p-values with both df multiplied by
#' epsilon (Huynh-Feldt epsilon capped at 1), as `anova.mlm(test =
#' "Spherical")`. The mixed ANOVA rows carry the same columns for the within
#' term and the interaction (epsilon from the pooled within-group covariance;
#' `NA` for the between term).
#' @inherit testrepeated details
#' @export
testrepeatedlong <- function(data, outcome, within, id, between = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE,
                             xlab, ylab, legendtitle,
                             posthoc = c("auto", "t", "wilcoxon", "nemenyi", "conover", "none"),
                             p.adjust.method = "BH") {
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  within_nm <- rlang::as_name(rlang::ensym(within))
  id_nm <- rlang::as_name(rlang::ensym(id))
  between_nm <- if (is.null(substitute(between))) NULL else rlang::as_name(rlang::ensym(between))
  out <- repeated_core(data, outcome_nm, within_nm, id_nm, outcome_label = outcome_nm, display_by_within = FALSE,
                       alpha = alpha, plot = plot, na.rm = na.rm, call = match.call(), between_nm = between_nm,
                       posthoc = match.arg(posthoc), p.adjust.method = p.adjust.method)
  out$plot <- relabel(out$plot, if (!missing(xlab)) list(x = xlab), if (!missing(ylab)) list(y = ylab),
                      if (!missing(legendtitle)) list(colour = legendtitle, shape = legendtitle))
  out
}

repeated_core <- function(data, outcome_nm, within_nm, id_nm, outcome_label = NULL, display_by_within = FALSE, alpha = 0.05, plot = TRUE, na.rm = TRUE, call = NULL, between_nm = NULL,
                          posthoc = "auto", p.adjust.method = "BH") {
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
  cc <- stats::complete.cases(wide[, measure_nms, drop = FALSE])
  complete <- wide[cc, measure_nms, drop = FALSE]
  sphericity <- check_sphericity(as.matrix(complete), alpha)
  friedman <- stats::friedman.test(as.matrix(complete))
  eps <- sphericity_eps(as.matrix(complete))
  anova <- add_eps(repeated_anova_test(df, outcome_nm, within_nm, id_nm), eps)
  if (ncol(complete) > 2) {
    sphericity$details <- paste0("Greenhouse-Geisser epsilon = ", format_stat(eps[["gg"]]),
                                 "; Huynh-Feldt epsilon = ", format_stat(eps[["hf"]]))
  }
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
  if (posthoc == "auto") posthoc <- if (recommendation == "Friedman test") "wilcoxon" else "t"
  posthoc_tests <- switch(posthoc, none = NULL,
    t = , wilcoxon = paired_posthoc_numeric(complete, c(t = "t", wilcoxon = "wilcox")[[posthoc]], p.adjust.method),
    friedman_pairs(complete, posthoc, p.adjust.method))
  alt <- list(anova = add_null_hypothesis(anova, h0), friedman = add_null_hypothesis(safe_tidy_htest(friedman, "Friedman test"), h0))
  if (!is.null(between_nm)) {
    bw <- df[[between_nm]][match(wide[[id_nm]][cc], df[[id_nm]])]
    alt$mixed_anova <- mixed_anova_test(df, outcome_nm, within_nm, id_nm, between_nm)
    m <- alt$mixed_anova
    eb <- sphericity_eps(as.matrix(complete), bw)
    alt$mixed_anova <- add_eps(m, eb, df2 = m$den.df, rows = m$term != between_nm)
  }
  plt <- if (plot) {
    sub <- if (is.null(between_nm)) plot_subtitle(recommendation, primary[1, ]) else termp("Mixed ANOVA", alt$mixed_anova)
    lab <- if (display_by_within) c("Measurement", "Value") else c(within_nm, outcome_nm)
    profile_plot(df, outcome_nm, within_nm, between_nm, id_nm, subtitle = sub, xlab = lab[1], ylab = lab[2],
                 caption = "Faint lines: subjects; bold: mean and 95% CI.")
  } else NULL
  out <- new_testflow("repeated", "repeated numeric measurements", outcome_label, within_nm, id_nm, df, descriptives_numeric(df, outcome_nm, within_nm), assumption_checks(check_independence_note("Repeated measurements from the same subjects are assumed by design."), normality, sphericity), list(test = recommendation, rationale = recommend_rationale(recommendation, normality, alpha = alpha)), primary, alt, posthoc = posthoc_tests, effect_size = effect, plot = plt, call = call, subclass = "repeated")
  out$interpretation <- make_report(out, alpha)
  if (recommendation == "Repeated-measures ANOVA" && ncol(complete) > 2) {
    out$interpretation <- paste0(
      out$interpretation, " Sphericity corrections: Greenhouse-Geisser epsilon = ", format_stat(eps[["gg"]]),
      ", p = ", format_p(anova$p.gg), "; Huynh-Feldt epsilon = ", format_stat(eps[["hf"]]),
      ", p = ", format_p(anova$p.hf), "."
    )
  }
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
