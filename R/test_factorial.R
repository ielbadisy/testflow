#' Run a factorial ANOVA workflow
#' @param formula A formula such as `outcome ~ factor1 * factor2`, or a data frame when using data-first style.
#' @param data A data frame, or the outcome column when using data-first style.
#' @param factors Factor columns selected with tidyselect syntax. Optional when using formula style.
#' @param alpha Significance level.
#' @param type Sums-of-squares type: `1` (sequential, base `aov()`), `2`, or
#'   `3` (via `car::Anova()`, required for `type = 2`/`3`). For unbalanced
#'   designs these can give materially different p-values; `type = 2` is a
#'   reasonable default when there is no strong prior reason to test factors
#'   in a particular order. `type = 3` requires sum-to-zero contrasts, which
#'   this function sets automatically.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_factorial`. The object is a
#' list containing the cleaned data, descriptive statistics, residual and
#' variance assumption checks, recommended factorial ANOVA, primary results,
#' ANOVA table, effect size, optional `ggplot`, original call, and report text.
#' `primary_test` has one row per model term (main effects and interactions)
#' with columns `term`, `statistic` (F), `parameter` and `num.df` (term df),
#' `den.df` (residual df), `p.value`, the term's eta squared, and its null
#' hypothesis. Levene's test for variance homogeneity is computed across all
#' design cells (every combination of factor levels).
#' @references
#' Fisher, R. A. (1925). \emph{Statistical Methods for Research Workers}.
#' Oliver and Boyd.
#'
#' Cohen, J. (1988). \emph{Statistical Power Analysis for the Behavioral
#' Sciences} (2nd ed.). Lawrence Erlbaum.
#' @export
testfactorial <- function(formula, data, factors = NULL, alpha = 0.05, type = 2, plot = TRUE, na.rm = TRUE) {
  first_expr <- tf_formula_arg(substitute(formula), parent.frame())
  second_expr <- tf_formula_arg(substitute(data), parent.frame())
  data_obj <- resolve_data_first_or_formula(formula, data)
  if (inherits(first_expr, "formula") || (is.call(first_expr) && identical(first_expr[[1]], as.name("~")))) {
    outcome_nm <- all.vars(first_expr[[2]])[1]
    factor_nms <- all.vars(first_expr[[3]])
    formula_obj <- stats::as.formula(paste(deparse(first_expr), collapse = " "))
  } else if (inherits(second_expr, "formula") || (is.call(second_expr) && identical(second_expr[[1]], as.name("~")))) {
    outcome_nm <- all.vars(second_expr[[2]])[1]
    factor_nms <- all.vars(second_expr[[3]])
    formula_obj <- stats::as.formula(paste(deparse(second_expr), collapse = " "))
  } else {
    outcome_nm <- rlang::as_name(second_expr)
    factor_nms <- tidyselect_names(data_obj, {{ factors }})
    formula_obj <- stats::as.formula(paste(tf_backtick(outcome_nm), "~", paste(tf_backtick(factor_nms), collapse = " * ")))
  }
  warn_if_two_groups_for_factorial(factor_nms)
  stopifnot(type %in% c(1, 2, 3))
  df <- drop_missing(data_obj, c(outcome_nm, factor_nms), na.rm = na.rm)
  df[factor_nms] <- lapply(df[factor_nms], as.factor)
  fit <- stats::aov(formula_obj, data = df)
  ss_tab <- if (identical(type, 1)) {
    stats::anova(fit)
  } else {
    if (identical(type, 3)) {
      # Type III SS are only interpretable with sum-to-zero contrasts;
      # under the default treatment contrasts they test a different (usually
      # meaningless) hypothesis about a reference-level intercept.
      old_contrasts <- options(contrasts = c("contr.sum", "contr.poly"))
      on.exit(options(old_contrasts), add = TRUE)
      lm_fit <- stats::lm(formula_obj, data = df)
    } else {
      lm_fit <- stats::lm(formula_obj, data = df)
    }
    car::Anova(lm_fit, type = type)
  }
  tab <- broom::tidy(ss_tab)
  tab <- tab[tab$term != "(Intercept)", , drop = FALSE]
  residual_df <- .bt(.resid = stats::residuals(fit))
  normality <- .residual_normality(residual_df$.resid, alpha, "Residuals appear approximately normal.", "Residuals deviate from normality.")
  # Levene's test across all design cells (every factor-level combination),
  # not only the levels of the first factor.
  cell_df <- df
  cell_df$.cell <- interaction(df[factor_nms], drop = TRUE, sep = ":")
  levene_test <- check_variance_homogeneity(cell_df, outcome_nm, ".cell", alpha)
  levene <- assumption_check(
    "Variance homogeneity",
    levene_test$status[1],
    ifelse(levene_test$status[1] == "acceptable", "Variance homogeneity looks reasonable.", "Variance homogeneity may be violated."),
    method = levene_test$method[1],
    statistic = levene_test$statistic[1],
    p_value = levene_test$p[1],
    details = paste0("Across ", nlevels(cell_df$.cell), " cells of ", paste(factor_nms, collapse = " x "), "; Df1=", levene_test$df1[1], "; Df2=", levene_test$df2[1])
  )
  type_label <- c(`1` = "Type I (sequential)", `2` = "Type II", `3` = "Type III")[[as.character(type)]]
  balanced <- assumption_check("Balanced design", "not required", ifelse(length(unique(table(df[factor_nms]))) > 1, paste0("Cell sizes are unbalanced; ", type_label, " sums of squares are used so the terms are not order-dependent.", if (identical(type, 1)) " Type I sums of squares depend on the order factors are listed in the formula for unbalanced designs; consider `type = 2` or `type = 3`." else ""), "Cell sizes are balanced; Type I, II, and III sums of squares agree."))
  effect <- eta_squared_aov(fit, tab = ss_tab)
  keep_rows <- tab[!is.na(tab$term) & tab$term != "Residuals", , drop = FALSE]
  resid_df <- tab$df[!is.na(tab$term) & tab$term == "Residuals"][1]
  primary <- .bt(
    method = paste0("Factorial ANOVA (", type_label, ")"),
    term = keep_rows$term,
    statistic = keep_rows$statistic,
    parameter = keep_rows$df,
    p.value = keep_rows$p.value,
    num.df = keep_rows$df,
    den.df = rep(resid_df, nrow(keep_rows))
  )
  plt <- if (plot && length(factor_nms) >= 2) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[factor_nms[1]]], y = .data[[outcome_nm]], color = .data[[factor_nms[2]]], group = .data[[factor_nms[2]]])) +
      ggplot2::stat_summary(fun = mean, geom = "line") +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 2) +
      ggplot2::labs(title = "Factorial ANOVA workflow", subtitle = plot_subtitle("Factorial ANOVA", list(p.value = min(primary$p.value, na.rm = TRUE))), x = factor_nms[1], y = outcome_nm, color = factor_nms[2]) +
      ggplot2::theme_minimal()
  } else if (plot) make_plot("groups", df, outcome_nm, factor_nms[1], "Factorial ANOVA", list(p.value = min(primary$p.value, na.rm = TRUE)), effect) else NULL
  h0 <- h0_factorial_terms(outcome_nm, primary$term)
  out <- new_testflow("factorial", "factorial design", outcome_nm, paste(factor_nms, collapse = ", "), data = df, descriptives = descriptives_numeric(df, outcome_nm, factor_nms[1]), assumptions = assumption_checks(check_independence_note(), normality, levene, balanced), recommended = list(test = "Factorial ANOVA", rationale = "Primary workflow for factorial numeric outcomes."), primary_test = add_null_hypothesis(primary, h0), alternative_tests = list(anova_table = tab), effect_size = effect, plot = plt, call = match.call(), subclass = "factorial")
  # One row per ANOVA term: attach each term's own eta squared.
  idx <- match(out$primary_test$term, effect$term)
  out$primary_test$effect_size <- effect$estimate[idx]
  out$primary_test$effect_size_name <- effect$name[idx]
  out$primary_test$effect_size_magnitude <- effect$magnitude[idx]
  out$interpretation <- make_report(out, alpha)
  out
}

h0_factorial_terms <- function(outcome, terms) {
  vapply(terms, function(tm) {
    parts <- strsplit(tm, ":", fixed = TRUE)[[1]]
    if (length(parts) > 1) {
      paste0("H0: there is no ", paste(parts, collapse = " x "), " interaction on the mean of ", outcome, ".")
    } else {
      paste0("H0: the population mean of ", outcome, " is equal across levels of ", tm, " (main effect).")
    }
  }, character(1), USE.NAMES = FALSE)
}

# Per-term summary used by print() and report() for factorial ANOVA.
factorial_terms_text <- function(x, alpha = NULL) {
  t <- x$primary_test
  vapply(seq_len(nrow(t)), function(i) {
    paste0(
      t$term[i], ": F = ", format_stat(t$statistic[i]),
      ", df = ", format_stat(t$num.df[i]), ", ", format_stat(t$den.df[i]),
      ", p = ", format_p(t$p.value[i]),
      if (!is.na(t$effect_size[i])) paste0(", ", t$effect_size_name[i], " = ", format_stat(t$effect_size[i])) else "",
      if (!is.null(alpha)) ifelse(t$p.value[i] < alpha, " (significant)", " (not significant)") else ""
    )
  }, character(1))
}
