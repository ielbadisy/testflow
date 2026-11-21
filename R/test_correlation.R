#' Test correlation between two numeric variables
#' @param formula A formula such as `y ~ x`, or a data frame when using pipe/data-first style.
#' @param data A data frame, or a first numeric column when using data-first style.
#' @param y Second numeric column. Optional when using formula style.
#' @param method Correlation method or `"auto"`. Only `"auto"` runs the
#'   Shapiro-Wilk pre-test (to choose Pearson or Spearman); a given method
#'   reports normality as "not checked". Shapiro-Wilk is never run for
#'   n > 5000 (normality is then "not checked"; use a Q-Q plot).
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param covariates Optional character vector of columns to partial out.
#'   When given, the primary result is the partial correlation of `y` and
#'   `x` given the covariates: the Pearson correlation of the residuals of
#'   `y` and `x` regressed on the covariates (ranks of all variables for
#'   Spearman). Its t test has \eqn{n - 2 - k} df for \eqn{k} covariate
#'   columns and the Fisher-z interval uses SE \eqn{1/\sqrt{n - 3 - k}}.
#'   The marginal result moves to `$alternative_tests$marginal`. Kendall is
#'   not supported. Default `NULL` keeps the marginal analysis.
#' @return A `testflow` object with class `testflow_correlation`. The object is
#' a list containing the cleaned complete-case data, numeric descriptives,
#' assumption checks, recommended correlation method, primary correlation test
#' with null hypothesis, Pearson/Spearman/Kendall results, a correlation table,
#' effect size, optional `ggplot`, original call, and report text.
#' @references
#' Pearson, K. (1895). Notes on regression and inheritance in the case of two
#' parents. \emph{Proceedings of the Royal Society of London}, 58, 240-242.
#'
#' Spearman, C. (1904). The proof and measurement of association between two
#' things. \emph{The American Journal of Psychology}, 15(1), 72-101.
#'
#' Kendall, M. G. (1938). A new measure of rank correlation. \emph{Biometrika},
#' 30(1/2), 81-93.
#' @export
testcorrelation <- function(formula, data, y = NULL, method = c("auto", "pearson", "spearman", "kendall"), alpha = 0.05, plot = TRUE, na.rm = TRUE,
                            covariates = NULL) {
  method <- match.arg(method)
  vars <- resolve_formula_pair(tf_formula_arg(substitute(formula), parent.frame()), tf_formula_arg(substitute(data), parent.frame()), substitute(y), missing(y), rhs_label = "x")
  y_nm <- vars$outcome
  x_nm <- vars$group
  data_obj <- resolve_data_first_or_formula(formula, data)
  df <- drop_missing(data_obj, c(x_nm, y_nm, covariates), na.rm = na.rm)
  warn_if_small_correlation_n(df, x_nm, y_nm)
  # The normality pre-test only drives the automatic choice; skip it when the method is given.
  normality <- check_normality(df, c(x_nm, y_nm), alpha = alpha, skip = if (method != "auto")
    paste0("Not pre-tested: method = \"", method, "\" was given; inspect a normal Q-Q plot if relevant."))
  outlier_flags <- iqr_outlier_flags(df, c(x_nm, y_nm))
  outliers <- iqr_outlier_assumption(outlier_flags, c(x_nm, y_nm))
  chosen <- if (method == "auto") {
    if (all(normality$status == "acceptable") && !any(outlier_flags$is_outlier, na.rm = TRUE)) "pearson" else "spearman"
  } else method
  linearity <- check_linearity(df[[x_nm]], df[[y_nm]], alpha)
  monotonicity <- check_monotonicity(df[[x_nm]], df[[y_nm]], alpha = alpha)
  pearson <- stats::cor.test(df[[x_nm]], df[[y_nm]], method = "pearson")
  spearman <- suppressWarnings(stats::cor.test(df[[x_nm]], df[[y_nm]], method = "spearman", exact = FALSE))
  kendall <- suppressWarnings(stats::cor.test(df[[x_nm]], df[[y_nm]], method = "kendall", exact = FALSE))
  primary <- switch(chosen, pearson = pearson, spearman = spearman, kendall = kendall)
  chosen_label <- title_case_method(paste(chosen, "correlation"))
  effect <- .bt(name = paste0(chosen_label, " r"), estimate = unname(primary$estimate), magnitude = magnitude_cramers_v(abs(unname(primary$estimate))))
  plt <- if (plot) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_nm]], y = .data[[y_nm]])) +
      ggplot2::geom_point(alpha = 0.6, colour = TFCOL[["pt"]]) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = TFCOL[["a1"]],
                           fill = TFCOL[["line"]]) +
      ggplot2::labs(subtitle = paste0(effstr(.bt(name = c(pearson = "Pearson r", spearman = "Spearman rho", kendall = "Kendall tau")[[chosen]],
                                                 estimate = effect$estimate[1], conf.low = primary$conf.int[1] %||% NA,
                                                 conf.high = primary$conf.int[2] %||% NA)), ", ", pstr(primary$p.value)),
                    caption = "Line: least-squares fit with 95% band.", x = x_nm, y = y_nm) +
      tftheme()
  } else NULL
  h0 <- h0_no_correlation(x_nm, y_nm)
  corr_assumptions <- switch(chosen, pearson = assumption_checks(linearity, normality, outliers), spearman = assumption_checks(monotonicity, outliers, assumption_check("Normality", "not required", "Normality is not required for Spearman correlation.")), kendall = assumption_checks(monotonicity, outliers, assumption_check("Normality", "not required", "Normality is not required for Kendall correlation.")))
  out <- new_testflow("correlation", "two numeric variables", y_nm, x_nm, data = df, descriptives = descriptives_numeric(df, c(x_nm, y_nm)), assumptions = corr_assumptions, recommended = list(test = chosen_label, rationale = if (method == "auto") recommend_rationale(chosen_label, normality, alpha = alpha, extra = if (any(outlier_flags$is_outlier, na.rm = TRUE)) "IQR outliers were flagged, which also favours Spearman" else "no IQR outliers were flagged") else paste0(chosen_label, " was requested by the user (method = \"", method, "\").")), primary_test = add_null_hypothesis(safe_tidy_htest(primary, chosen_label), h0), alternative_tests = list(correlation_table = .bt(method = c("Pearson", "Spearman", "Kendall"), statistic = c(unname(pearson$estimate), unname(spearman$estimate), unname(kendall$estimate)), p.value = c(pearson$p.value, spearman$p.value, kendall$p.value)), pearson = add_null_hypothesis(safe_tidy_htest(pearson, "Pearson correlation"), h0), spearman = add_null_hypothesis(safe_tidy_htest(spearman, "Spearman correlation"), h0), kendall = add_null_hypothesis(safe_tidy_htest(kendall, "Kendall correlation"), h0)), effect_size = effect, plot = plt, call = match.call(), subclass = "correlation")
  if (length(covariates)) {
    if (chosen == "kendall") stop("Partial correlation supports Pearson or Spearman only.", call. = FALSE)
    pc <- partial_cor(df[[x_nm]], df[[y_nm]], df[covariates], chosen, alpha)
    lab <- paste0("Partial ", chosen_label, " (given ", paste(covariates, collapse = ", "), ")")
    h0p <- paste0("H0: ", x_nm, " and ", y_nm, " are uncorrelated given ", paste(covariates, collapse = ", "), ".")
    out$alternative_tests <- c(list(marginal = out$primary_test), out$alternative_tests)
    out$primary_test <- add_null_hypothesis(safe_tidy_htest(pc, lab), h0p)
    out$primary_test$conf.level <- 1 - alpha
    out$recommended$test <- lab
    out$effect_size <- .bt(name = paste0(lab, " r"), estimate = unname(pc$estimate),
                           magnitude = magnitude_cramers_v(abs(unname(pc$estimate))))
  }
  out$interpretation <- make_report(out, alpha)
  out
}

iqr_outlier_assumption <- function(data, vars) {
  n <- sum(data$is_outlier, na.rm = TRUE)
  assumption_check(
    "Extreme outliers",
    ifelse(n == 0, "acceptable", "warning"),
    paste0(n, " potential outlier(s) flagged by IQR."),
    details = paste0("IQR rule applied to ", paste(vars, collapse = ", "))
  )
}

iqr_outlier_flags <- function(data, vars) {
  .map_dfr(vars, function(v) {
    x <- data[[v]]
    x <- x[!is.na(x)]
    if (length(x) < 4) {
      return(.bt(variable = v, is_outlier = NA))
    }
    q1 <- stats::quantile(x, 0.25, names = FALSE)
    q3 <- stats::quantile(x, 0.75, names = FALSE)
    i <- q3 - q1
    .bt(variable = v, is_outlier = x < q1 - 1.5 * i | x > q3 + 1.5 * i)
  })
}

#' Test a correlation matrix
#' @param data A data frame.
#' @param vars Numeric columns.
#' @param method Correlation method.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param covariates Optional character vector of columns to partial out of
#'   every pair (see [testcorrelation()]); each pair then uses its rows
#'   complete on the pair and the covariates. Pearson or Spearman only.
#' @return A `testflow` object with class `testflow_correlation_matrix`. The
#' object is a list containing the cleaned data, numeric descriptives,
#' screening assumptions, selected correlation-matrix method, pairwise
#' correlation matrix, pairwise p-value table, maximum absolute correlation as
#' an effect-size summary, optional heatmap `ggplot`, original call, and report
#' text. The pairwise results (one row per variable pair: `var1`, `var2`,
#' `estimate`, raw `p`, pairwise-complete `n`, and Benjamini-Hochberg adjusted
#' `p.adj`) are in `$posthoc`, and also in
#' `$alternative_tests$correlation_table`; the matrix is in
#' `$alternative_tests$correlation_matrix`. `$primary_test` is only a summary:
#' its `p.value` is the smallest unadjusted pairwise p-value and its
#' `statistic` is `NA`. Use `$posthoc$p.adj` for inference on single pairs.
#' @export
testcorrelationmatrix <- function(data, vars, method = c("spearman", "pearson", "kendall"), alpha = 0.05, plot = TRUE,
                                  na.rm = TRUE, covariates = NULL) {
  method <- match.arg(method)
  vars <- tidyselect_names(data, {{ vars }})
  warn_if_screening_workflow("correlation_matrix")
  df <- drop_missing(data, vars, na.rm = na.rm)
  # Each pair's correlation/test uses its own pairwise-complete rows, taken
  # from the un-listwise-filtered columns; listwise-dropping across all
  # `vars` first (as `df` does, for the stored/descriptive data) would
  # needlessly discard rows for pairs that don't involve the variable
  # carrying the missing value, and would make `use = "pairwise.complete.obs"`
  # below a no-op.
  full_vars <- as.data.frame(data)[, vars, drop = FALSE]
  mat <- stats::cor(full_vars, method = method, use = "pairwise.complete.obs")
  pairs <- utils::combn(vars, 2, simplify = FALSE)
  if (length(covariates) && method == "kendall") stop("Partial correlation supports Pearson or Spearman only.", call. = FALSE)
  cov_df <- as.data.frame(data)[, covariates, drop = FALSE]
  long <- .map_dfr(pairs, function(z) {
    if (length(covariates)) {
      ok <- stats::complete.cases(full_vars[z], cov_df)
      test <- partial_cor(full_vars[[z[1]]][ok], full_vars[[z[2]]][ok], cov_df[ok, , drop = FALSE], method, alpha)
      return(.bt(var1 = z[1], var2 = z[2], estimate = unname(test$estimate), p = test$p.value, n = sum(ok)))
    }
    test <- suppressWarnings(stats::cor.test(full_vars[[z[1]]], full_vars[[z[2]]], method = method, exact = FALSE))
    .bt(var1 = z[1], var2 = z[2], estimate = unname(test$estimate), p = test$p.value, n = sum(stats::complete.cases(full_vars[[z[1]]], full_vars[[z[2]]])))
  })
  if (length(covariates)) {
    mat[] <- 1
    for (i in seq_len(nrow(long))) mat[long$var1[i], long$var2[i]] <- mat[long$var2[i], long$var1[i]] <- long$estimate[i]
  }
  long$p.adj <- stats::p.adjust(long$p, method = "BH")
  long$p.adjust.method <- "BH"
  n_sig_adj <- sum(long$p.adj < alpha, na.rm = TRUE)
  assumptions <- assumption_checks(
    assumption_check("Variable types", "acceptable", "All selected variables should be numeric or ordinal."),
    assumption_check("Pairwise missing data", "acceptable", "Pairwise complete observations are used for each correlation."),
    assumption_check(
      "Multiple testing correction", "acceptable",
      paste0(n_sig_adj, " of ", nrow(long), " pairwise correlation(s) remain significant after Benjamini-Hochberg adjustment (see the `p.adj` column)."),
      method = "Benjamini-Hochberg (BH)"
    )
  )
  plt <- if (plot) {
    hm <- as.data.frame(as.table(mat))
    hm$label <- sprintf("%.2f", hm$Freq)
    hm$text_color <- ifelse(abs(hm$Freq) >= 0.65, "light", "dark")
    ggplot2::ggplot(hm, ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$Freq)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.6) +
      ggplot2::geom_text(ggplot2::aes(label = .data$label, color = .data$text_color), size = 4) +
      ggplot2::scale_fill_gradient2(limits = c(-1, 1)) +
      ggplot2::scale_color_manual(values = c(dark = "#1F2933", light = "white"), guide = "none") +
      ggplot2::coord_equal() +
      ggplot2::labs(title = "Correlation matrix", x = NULL, y = NULL, fill = "r") +
      ggplot2::theme_minimal()
  } else NULL
  method_label <- title_case_method(paste(method, "correlation matrix"))
  primary <- .bt(method = method_label, statistic = NA_real_, parameter = NA_real_, p.value = min(long$p, na.rm = TRUE))
  out <- new_testflow("correlation_matrix", "correlation matrix", paste(vars, collapse = ", "), data = df, descriptives = descriptives_numeric(df, vars), assumptions = assumptions, recommended = list(test = method_label), primary_test = primary, alternative_tests = list(correlation_matrix = mat, correlation_table = long), posthoc = long, effect_size = .bt(name = "Maximum absolute r", estimate = max(abs(long$estimate), na.rm = TRUE), magnitude = magnitude_cramers_v(max(abs(long$estimate), na.rm = TRUE))), plot = plt, call = match.call(), subclass = "correlation_matrix")
  out$interpretation <- make_report(out, alpha)
  out
}

tidyselect_names <- function(data, expr) {
  if (rlang::quo_is_symbol(rlang::enquo(expr))) {
    nm <- rlang::as_name(rlang::ensym(expr))
    if (nm %in% names(data)) return(nm)
  }
  select_names(data, rlang::enquo(expr))
}

# Base-R column selection for the tidyselect forms used with testflow: bare
# names, strings, c(), a:b ranges, -x exclusions, character vectors from the
# calling environment, and starts_with(), ends_with(), contains(), matches(),
# everything(), all_of() and any_of().
select_names <- function(data, quo) {
  nms <- names(data)
  expr <- rlang::quo_get_expr(quo)
  env <- rlang::quo_get_env(quo)
  if (is.null(expr)) return(character())
  mask <- new.env(parent = env)
  for (i in seq_along(nms)) assign(nms[i], i, envir = mask)
  pos <- function(x) which(x)
  mask$starts_with <- function(match, ignore.case = TRUE) {
    pos(startsWith(if (ignore.case) tolower(nms) else nms, if (ignore.case) tolower(match) else match))
  }
  mask$ends_with <- function(match, ignore.case = TRUE) {
    pos(endsWith(if (ignore.case) tolower(nms) else nms, if (ignore.case) tolower(match) else match))
  }
  mask$contains <- function(match, ignore.case = TRUE) {
    pos(grepl(match, nms, fixed = !ignore.case, ignore.case = ignore.case))
  }
  mask$matches <- function(match, ignore.case = TRUE, perl = FALSE) {
    pos(grepl(match, nms, ignore.case = ignore.case, perl = perl))
  }
  mask$everything <- function() seq_along(nms)
  mask$all_of <- function(x) {
    miss <- setdiff(x, nms)
    if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
    match(x, nms)
  }
  mask$any_of <- function(x) match(intersect(x, nms), nms)
  val <- eval(expr, mask)
  if (is.character(val)) {
    miss <- setdiff(val, nms)
    if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
    return(unique(val))
  }
  if (!is.numeric(val)) stop("Column selection must give column names or positions.", call. = FALSE)
  val <- as.integer(val)
  if (any(abs(val) > length(nms) | val == 0L)) stop("Column selection is out of range.", call. = FALSE)
  keep <- unique(val[val > 0L])
  drop <- -val[val < 0L]
  if (!length(keep) && length(drop)) keep <- seq_along(nms)
  nms[setdiff(keep, drop)]
}

# Partial correlation as an htest: correlation of residuals on the covariates (on ranks for Spearman).
partial_cor <- function(x, y, z, method = "pearson", alpha = 0.05) {
  rk <- function(v) if (method == "spearman") rank(v) else v
  Z <- stats::model.matrix(~ ., data.frame(lapply(z, function(v) if (is.numeric(v)) rk(v) else v)))
  r <- stats::cor(stats::lm.fit(Z, rk(x))$residuals, stats::lm.fit(Z, rk(y))$residuals)
  n <- length(x); k <- ncol(Z) - 1; df <- n - 2 - k; tt <- r * sqrt(df / (1 - r^2))
  ci <- tanh(atanh(r) + c(-1, 1) * stats::qnorm(1 - alpha / 2) / sqrt(n - 3 - k))
  structure(list(statistic = c(t = tt), parameter = c(df = df), p.value = 2 * stats::pt(-abs(tt), df),
                 estimate = c(cor = r), conf.int = structure(ci, conf.level = 1 - alpha), alternative = "two.sided",
                 method = paste0("Partial ", method, " correlation")), class = "htest")
}
