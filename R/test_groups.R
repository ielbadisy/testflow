#' Compare a numeric outcome across more than two groups
#' @param formula A formula such as `outcome ~ group`, or a data frame when using pipe/data-first style.
#' @param data A data frame, or an outcome column when using data-first style.
#' @param group Grouping column. Optional when using formula style.
#' @param alpha Significance level.
#' @param posthoc Post-hoc method. `"auto"` (default, same as `TRUE`) follows
#'   the recommended omnibus test: Tukey HSD after one-way ANOVA, pairwise
#'   Welch t tests after Welch ANOVA, pairwise Wilcoxon rank-sum tests after
#'   Kruskal-Wallis. `"tukey"`, `"gameshowell"` (Games-Howell, unequal
#'   variances), `"welch"`, `"wilcoxon"`, `"dunn"` (Dunn's test on pooled
#'   ranks, see [testdunn()]), `"conover"` (Conover-Iman test on pooled ranks)
#'   or `"dscf"` (Dwass-Steel-Critchlow-Fligner test) fixes the method;
#'   `"none"` (or `FALSE`) skips it. See Details for the rank tests.
#' @param p.adjust.method Multiplicity adjustment passed to
#'   [stats::p.adjust()] for the `"welch"`, `"wilcoxon"`, `"dunn"` and
#'   `"conover"` post hoc tests (Tukey, Games-Howell and DSCF carry their own
#'   studentized-range adjustment).
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_groups`. The object is a
#' list containing the cleaned data, grouped descriptive statistics, assumption
#' checks, recommended omnibus test, primary test result with null hypothesis,
#' alternative omnibus results, post-hoc comparisons when requested, effect
#' size, optional `ggplot`, original call, and report text.
#'
#' `$posthoc` is a list with `method`, `result` (a `TukeyHSD` or
#' `pairwise.htest` object), `p.adjust.method` and `table`, a tidy table with
#' one row per pair: `group1`, `group2`, `estimate` (group2 minus group1: mean
#' difference, mean rank difference for Dunn and Conover-Iman, Hodges-Lehmann shift
#' for DSCF), `statistic`, `conf.low`,
#' `conf.high`, `p`, `p.adj` and `p.adjust.method` (`NA` where a method does
#' not give the quantity).
#'
#' After a Kruskal-Wallis test, `$effect_size` has two rows:
#' `eta^2_H = (H - k + 1) / (n - k)` and `epsilon^2_H = H / (n - 1)`.
#' @details
#' Rank post hoc tests. Let \eqn{r} be the ranks of all \eqn{N} values (ties
#' get mid-ranks), \eqn{\bar R_j} and \eqn{n_j} the mean rank and size of
#' group \eqn{j}, and \eqn{k} the number of groups. Unequal group sizes enter
#' through \eqn{n_i} and \eqn{n_j}.
#'
#' Conover-Iman (1979): with the rank variance
#' \eqn{S^2 = (\sum r^2 - N(N+1)^2/4)/(N-1)} and \eqn{H} the tie-corrected
#' Kruskal-Wallis statistic,
#' \deqn{t_{ij} = \frac{\bar R_j - \bar R_i}{\sqrt{S^2 \frac{N-1-H}{N-k}
#'   \left(\frac{1}{n_i} + \frac{1}{n_j}\right)}},}
#' referred to a t distribution with \eqn{N-k} df, then adjusted by
#' `p.adjust.method`. \eqn{S^2} uses the observed ranks, which corrects for
#' ties. `estimate` is \eqn{\bar R_j - \bar R_i}.
#'
#' Dwass-Steel-Critchlow-Fligner (Critchlow and Fligner 1991): each pair is
#' ranked on its own. With \eqn{U_j} the Mann-Whitney count of group \eqn{j}
#' over group \eqn{i}, \eqn{s = n_i + n_j} and \eqn{t} the tie sizes in the
#' pair,
#' \deqn{W_{ij} = \sqrt{2}\,\frac{U_j - n_i n_j / 2}{\sqrt{\frac{n_i n_j}{s(s-1)}
#'   \left(\frac{s^3 - s}{12} - \sum \frac{t^3 - t}{12}\right)}},}
#' with p from the studentized range with \eqn{k} means and infinite df
#' (single step, `p.adj = p`). `estimate` is the Hodges-Lehmann shift, the
#' median of all differences between group2 and group1 values.
#'
#' The Dunn, Conover-Iman, DSCF and Games-Howell tables match
#' `PMCMRplus::kwAllPairsDunnTest()`, `kwAllPairsConoverTest()`,
#' `dscfAllPairsTest()` and `gamesHowellTest()`. Games-Howell needs two
#' values per group (pairs with a singleton group are `NA`); pairwise Welch t
#' tests stop in that case. Use `plot(x, type = "posthoc")` to draw the table.
#' @references
#' Conover, W. J., and Iman, R. L. (1979). On multiple-comparisons procedures.
#' Technical Report LA-7677-MS, Los Alamos Scientific Laboratory.
#'
#' Critchlow, D. E., and Fligner, M. A. (1991). On distribution-free multiple
#' comparisons in the one-way analysis of variance. \emph{Communications in
#' Statistics - Theory and Methods}, 20(1), 127-139.
#'
#' Dunn, O. J. (1964). Multiple comparisons using rank sums.
#' \emph{Technometrics}, 6(3), 241-252.
#' @export
testgroups <- function(formula, data, group = NULL, alpha = 0.05,
                       posthoc = c("auto", "tukey", "gameshowell", "welch", "wilcoxon", "dunn", "conover", "dscf",
                                   "none"),
                       plot = TRUE, na.rm = TRUE, p.adjust.method = "BH") {
  posthoc <- if (isTRUE(posthoc)) "auto" else if (isFALSE(posthoc)) "none" else match.arg(posthoc)
  vars <- resolve_formula_pair(tf_formula_arg(substitute(formula), parent.frame()), tf_formula_arg(substitute(data), parent.frame()), substitute(group), missing(group))
  outcome_nm <- vars$outcome
  group_nm <- vars$group
  data_obj <- resolve_data_first_or_formula(formula, data)
  df <- drop_missing(data_obj, c(outcome_nm, group_nm), na.rm = na.rm)
  df[[group_nm]] <- droplevels(as.factor(df[[group_nm]]))
  if (length(unique(df[[group_nm]])) < 3) stop("`group` must contain at least three groups.", call. = FALSE)
  normality <- check_normality(df, outcome_nm, group_nm, alpha)
  levene <- check_variance_homogeneity(df, outcome_nm, group_nm, alpha)
  bartlett <- check_bartlett(df, outcome_nm, group_nm, alpha)
  independence <- check_independence_note()
  outliers <- check_outliers(df[[outcome_nm]])
  recommendation <- recommend_groups(normality, levene)
  formula <- tf_formula(outcome_nm, group_nm)
  aov_fit <- stats::aov(formula, data = df)
  aov_full <- broom::tidy(aov_fit)
  aov_test <- aov_full[1, , drop = FALSE]
  aov_test$num.df <- aov_full$df[1]
  aov_test$den.df <- aov_full$df[nrow(aov_full)]
  # Welch ANOVA needs a variance in every group; with a singleton group it is reported as NA.
  welch <- if (all(table(df[[group_nm]]) > 1)) stats::oneway.test(formula, data = df, var.equal = FALSE) else
    structure(list(statistic = c(F = NA_real_), parameter = c("num df" = NA_real_, "denom df" = NA_real_),
                   p.value = NA_real_, method = "One-way analysis of means (not assuming equal variances)",
                   data.name = deparse(formula)), class = "htest")
  kruskal <- stats::kruskal.test(formula, data = df)
  primary <- switch(recommendation, "One-way ANOVA" = list_obj_from_tidy(aov_test, "One-way ANOVA"), "Welch ANOVA" = welch, "Kruskal-Wallis test" = kruskal)
  h0 <- h0_mean_equal(outcome_nm, group_nm)
  primary_tidy <- if (inherits(primary, "htest")) safe_tidy_htest(primary, recommendation) else primary
  primary_tidy <- add_null_hypothesis(primary_tidy, h0)
  ph <- if (posthoc != "none") {
    posthoc_groups(df, outcome_nm, group_nm, recommendation, alpha, posthoc, p.adjust.method)
  }
  effect <- if (recommendation == "Kruskal-Wallis test") {
    h <- unname(kruskal$statistic); n <- nrow(df); k <- length(unique(df[[group_nm]]))
    es <- c((h - k + 1) / (n - k), h / (n - 1))
    .bt(name = c("Kruskal eta squared (H)", "Kruskal epsilon squared (H)"), estimate = es,
        magnitude = magnitude_eta2(es))
  } else eta_squared_aov(aov_fit)[1, , drop = FALSE]
  plt <- if (plot) make_plot("groups", df, outcome_nm, group_nm, recommendation, if (inherits(primary, "htest")) primary else list(p.value = primary_tidy$p.value[1]), effect) else NULL
  out <- new_testflow("groups", "more than two independent groups", outcome_nm, group_nm, data = df, descriptives = descriptives_numeric(df, outcome_nm, group_nm), assumptions = assumption_checks(independence, normality, levene, bartlett, outliers), recommended = list(test = recommendation, rationale = recommend_rationale(recommendation, normality, levene, alpha)), primary_test = primary_tidy, alternative_tests = list(anova = add_null_hypothesis(aov_test, h0), welch = add_null_hypothesis(safe_tidy_htest(welch, "Welch ANOVA"), h0), kruskal = add_null_hypothesis(safe_tidy_htest(kruskal, "Kruskal-Wallis test"), h0)), posthoc = ph, effect_size = effect, plot = plt, call = match.call(), subclass = "groups")
  out$interpretation <- make_report(out, alpha)
  out
}

list_obj_from_tidy <- function(x, method) {
  .bt(method = method, statistic = x$statistic[1], parameter = x$df[1], num.df = x[["num.df"]][1], den.df = x[["den.df"]][1], p.value = x$p.value[1])
}
