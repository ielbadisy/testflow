#' Test for trend across ordered categories
#'
#' Cochran-Armitage test for a binary outcome, or the linear-by-linear
#' (Mantel) test for an ordinal outcome, across the ordered levels of an
#' exposure, with the chi-square for departure from linear trend.
#'
#' @param formula `outcome ~ exposure`. The exposure levels are taken in
#'   factor order (or numeric order).
#' @param data A data frame.
#' @param scores Numeric scores for the exposure levels. Default `NULL` uses
#'   the values of a numeric exposure, otherwise `1, 2, ..., k`.
#' @param outcome_scores Numeric scores for the outcome levels (ordinal
#'   outcome). Default `NULL` uses numeric values or `1, 2, ...`.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_trend`. `$primary_test`
#'   holds the trend chi-square (`statistic`, `parameter` = 1, `p.value`) and
#'   `estimate`, the Pearson correlation \eqn{r} between the scores (its sign
#'   gives the direction). `$alternative_tests$trend_table` lists the trend,
#'   the departure from linear trend (`k - 2` df) and the overall Pearson
#'   chi-square; `$alternative_tests$by_level` gives counts and proportions
#'   (binary outcome) or mean outcome scores by exposure level.
#' @details
#' With \eqn{n} observations and \eqn{r} the correlation between the
#' exposure and outcome scores, the Cochran-Armitage statistic is
#' \eqn{n r^2} (as `stats::prop.trend.test()`) and the linear-by-linear
#' statistic is \eqn{(n - 1) r^2} (as `coin::lbl_test()`), both on 1 df.
#' @references
#' Armitage P. Tests for linear trends in proportions and frequencies.
#' *Biometrics*. 1955;11:375-386.
#'
#' Agresti A. *Categorical Data Analysis* (3rd ed.). Wiley; 2013. Section 3.4.
#' @examples
#' d <- data.frame(dose = rep(1:4, each = 40), y = rbinom(160, 1, rep(c(.1, .2, .3, .4), each = 40)))
#' testtrend(y ~ dose, d, plot = FALSE)
#' @export
testtrend <- function(formula, data, scores = NULL, outcome_scores = NULL, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  v <- all.vars(formula); y_nm <- v[1]; x_nm <- v[2]
  df <- drop_missing(data, c(y_nm, x_nm), na.rm = na.rm)
  lev_scores <- function(z, s) {
    f <- factor(z); lv <- levels(f)
    s <- s %||% if (is.numeric(z)) as.numeric(lv) else seq_along(lv)
    if (length(s) != length(lv)) stop("Scores must have one value per level.", call. = FALSE)
    list(f = f, s = s[as.integer(f)], lv = lv, sc = s)
  }
  ex <- lev_scores(df[[x_nm]], scores); oc <- lev_scores(df[[y_nm]], outcome_scores)
  k <- length(ex$lv); n <- nrow(df); binary <- length(oc$lv) == 2
  r <- stats::cor(ex$s, oc$s)
  trend <- (if (binary) n else n - 1) * r^2
  tab <- table(oc$f, ex$f)
  overall <- unname(suppressWarnings(stats::chisq.test(tab, correct = FALSE))$statistic)
  # Pearson X2 and n r^2 share the n denominator; LBL uses n - 1, so rescale the overall statistic to match.
  overall_s <- if (binary) overall else overall * (n - 1) / n
  method <- if (binary) "Cochran-Armitage test for trend" else "Linear-by-linear association test"
  pt <- stats::pchisq(trend, 1, lower.tail = FALSE)
  dep_df <- (k - 1) * (length(oc$lv) - 1) - 1
  trend_table <- .bt(component = c("Linear trend", "Departure from linear trend", "Overall association"),
                     statistic = c(trend, overall_s - trend, overall_s), parameter = c(1, dep_df, dep_df + 1),
                     p.value = c(pt, if (dep_df > 0) stats::pchisq(overall_s - trend, dep_df, lower.tail = FALSE) else NA,
                                 stats::pchisq(overall_s, dep_df + 1, lower.tail = FALSE)))
  by_level <- if (binary) {
    .bt(level = ex$lv, score = ex$sc, events = as.vector(tab[2, ]), n = as.vector(colSums(tab)),
        proportion = as.vector(tab[2, ] / colSums(tab)))
  } else .bt(level = ex$lv, score = ex$sc, n = as.vector(colSums(tab)),
             mean_outcome_score = as.vector(tapply(oc$s, ex$f, mean)))
  h0 <- paste0("H0: no linear trend in ", y_nm, " across the ordered levels of ", x_nm, ".")
  primary <- .bt(method = method, estimate = r, statistic = trend, parameter = 1, p.value = pt)
  plt <- if (plot) {
    pd <- data.frame(score = by_level$score, y = if (binary) by_level$proportion else by_level$mean_outcome_score)
    ggplot2::ggplot(pd, ggplot2::aes(x = .data$score, y = .data$y)) +
      ggplot2::geom_point(size = 3, colour = "#4C78A8") + ggplot2::geom_line(colour = "#4C78A8") +
      ggplot2::scale_x_continuous(breaks = by_level$score, labels = by_level$level) +
      ggplot2::labs(title = "Trend workflow", subtitle = paste0(method, ", p = ", format_p(pt)), x = x_nm,
                    y = if (binary) paste0("Proportion ", y_nm, " = ", oc$lv[2]) else paste0("Mean score of ", y_nm)) +
      ggplot2::theme_minimal()
  } else NULL
  out <- new_testflow("trend", "ordered exposure and categorical outcome", y_nm, x_nm, data = df,
                      descriptives = descriptives_categorical(df, c(y_nm, x_nm)),
                      assumptions = assumption_checks(check_independence_note(), assumption_check("Scores", "assumed",
                        paste0("Exposure scores: ", paste(ex$sc, collapse = ", "), "; the test targets a linear trend on these scores."))),
                      recommended = list(test = method, rationale = paste0(if (binary) "Binary" else "Ordinal",
                        " outcome across ordered exposure levels; one df targets the linear trend.")),
                      primary_test = add_null_hypothesis(primary, h0),
                      alternative_tests = list(trend_table = trend_table, by_level = by_level),
                      effect_size = .bt(name = "Correlation of scores (r)", estimate = r, magnitude = NA_character_),
                      plot = plt, call = match.call(), subclass = "trend")
  out$interpretation <- make_report(out, alpha)
  out
}

# Expected a and its variance in one 2 x 2 table under a common odds ratio (Breslow-Day).
bd_cell <- function(t, psi) {
  m1 <- sum(t[1, ]); n1 <- sum(t[, 1]); N <- sum(t)
  lo <- max(0, m1 + n1 - N); hi <- min(m1, n1)
  f <- function(a) a * (N - m1 - n1 + a) - psi * (m1 - a) * (n1 - a)
  a <- if (hi > lo) stats::uniroot(f, c(lo, hi), tol = 1e-12)$root else lo
  c(a = a, v = 1 / sum(1 / c(a, m1 - a, n1 - a, N - m1 - n1 + a)))
}

#' Stratified 2 x 2 analysis: Cochran-Mantel-Haenszel test and pooled odds ratio
#'
#' Tests the association between a binary exposure and a binary outcome
#' within strata, pools the stratum odds ratios (Mantel-Haenszel) and tests
#' their homogeneity (Breslow-Day, with Tarone's correction).
#'
#' @param formula `outcome ~ exposure | stratum`.
#' @param data A data frame.
#' @param event Outcome level counted as the event. Default `NULL` uses the
#'   second level.
#' @param correct Logical; continuity correction of the CMH statistic (as
#'   `stats::mantelhaen.test()`, default `TRUE`).
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_stratified`.
#'   `$primary_test` holds the CMH chi-square with the MH common odds ratio
#'   (`estimate`) and its Robins-Breslow-Greenland interval (`conf.low`,
#'   `conf.high`). The odds ratio compares the second exposure level with
#'   the first, for the `event` level. `$alternative_tests` has
#'   `stratum_or` (Woolf OR per stratum, 0.5 added in tables with a zero
#'   cell), `crude_or` (collapsed table), `breslow_day` and
#'   `breslow_day_tarone` (homogeneity, `K - 1` df).
#' @details
#' The Breslow-Day statistic is \eqn{\sum_k (a_k - E_k)^2 / V_k}, with
#' \eqn{E_k} and \eqn{V_k} the mean and asymptotic variance of the first
#' cell under the MH common odds ratio. Tarone's version subtracts
#' \eqn{(\sum a_k - \sum E_k)^2 / \sum V_k}. These match
#' `DescTools::BreslowDayTest()`.
#' @references
#' Mantel N, Haenszel W. Statistical aspects of the analysis of data from
#' retrospective studies of disease. *JNCI*. 1959;22:719-748.
#'
#' Breslow NE, Day NE. *Statistical Methods in Cancer Research, Vol. 1*.
#' IARC; 1980.
#'
#' Tarone RE. On heterogeneity tests based on efficient scores.
#' *Biometrika*. 1985;72:91-95.
#' @examples
#' d <- as.data.frame(UCBAdmissions)
#' d <- d[rep(seq_len(nrow(d)), d$Freq), 1:3]
#' teststratified(Admit ~ Gender | Dept, d, plot = FALSE)
#' @export
teststratified <- function(formula, data, event = NULL, correct = TRUE, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  v <- all.vars(formula); y_nm <- v[1]; x_nm <- v[2]; s_nm <- v[3]
  df <- drop_missing(data, v[1:3], na.rm = na.rm)
  y <- factor(df[[y_nm]]); x <- factor(df[[x_nm]]); s <- factor(df[[s_nm]])
  if (nlevels(y) != 2 || nlevels(x) != 2) stop("Outcome and exposure must both have two levels.", call. = FALSE)
  ev <- event %||% levels(y)[2]
  # Rows: exposure level 2 first; columns: event first, so a = exposed events and OR = odds(event | lvl2) / odds(event | lvl1).
  arr <- table(factor(x, rev(levels(x))), factor(y, c(ev, setdiff(levels(y), ev))), s)
  arr <- arr[, , apply(arr, 3, sum) > 0, drop = FALSE]
  K <- dim(arr)[3]
  mh <- stats::mantelhaen.test(arr, correct = correct, conf.level = 1 - alpha)
  psi <- unname(mh$estimate)
  bd <- vapply(seq_len(K), function(k) bd_cell(arr[, , k], psi), numeric(2))
  a <- arr[1, 1, ]
  bdx <- sum((a - bd["a", ])^2 / bd["v", ]); tar <- bdx - (sum(a) - sum(bd["a", ]))^2 / sum(bd["v", ])
  woolf <- function(t) {
    t <- t + if (any(t == 0)) 0.5 else 0; se <- sqrt(sum(1 / t)); z <- stats::qnorm(1 - alpha / 2)
    or <- t[1, 1] * t[2, 2] / (t[1, 2] * t[2, 1]); c(or, exp(log(or) + c(-1, 1) * z * se))
  }
  so <- vapply(seq_len(K), function(k) woolf(arr[, , k]), numeric(3))
  cr <- woolf(apply(arr, c(1, 2), sum))
  cmp <- paste0(levels(x)[2], " vs ", levels(x)[1])
  bdrow <- function(m, st) .bt(method = m, statistic = st, parameter = K - 1,
                               p.value = stats::pchisq(st, K - 1, lower.tail = FALSE))
  primary <- .bt(method = mh$method, estimate = psi, statistic = unname(mh$statistic), parameter = 1,
                 p.value = mh$p.value, conf.low = mh$conf.int[1], conf.high = mh$conf.int[2], conf.level = 1 - alpha)
  h0 <- paste0("H0: ", y_nm, " and ", x_nm, " are independent within every level of ", s_nm, " (common OR = 1).")
  stratum_or <- .bt(stratum = dimnames(arr)[[3]], n = as.vector(apply(arr, 3, sum)), events_exposed = as.vector(a),
                    odds_ratio = so[1, ], conf.low = so[2, ], conf.high = so[3, ])
  plt <- if (plot) {
    pd <- rbind(data.frame(stratum = stratum_or$stratum, or = so[1, ], lo = so[2, ], hi = so[3, ]),
                data.frame(stratum = "MH pooled", or = psi, lo = mh$conf.int[1], hi = mh$conf.int[2]))
    pd$stratum <- factor(pd$stratum, rev(pd$stratum))
    ggplot2::ggplot(pd, ggplot2::aes(x = .data$or, y = .data$stratum)) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
      ggplot2::geom_pointrange(ggplot2::aes(xmin = .data$lo, xmax = .data$hi), colour = "#4C78A8") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Stratified 2 x 2 workflow", subtitle = paste0("CMH p = ", format_p(mh$p.value),
                    "; Breslow-Day p = ", format_p(stats::pchisq(bdx, K - 1, lower.tail = FALSE))),
                    x = paste0("Odds ratio (", cmp, ")"), y = s_nm) +
      ggplot2::theme_minimal()
  } else NULL
  out <- new_testflow("stratified", "two binary variables within strata", y_nm, x_nm, data = df,
                      descriptives = descriptives_categorical(df, c(y_nm, x_nm, s_nm)),
                      assumptions = assumption_checks(check_independence_note(),
                        assumption_check("Homogeneity of odds ratios", if (stats::pchisq(bdx, K - 1, lower.tail = FALSE) < alpha)
                          "not acceptable" else "acceptable", "A common odds ratio across strata is assumed by the pooled estimate.",
                          method = "Breslow-Day", statistic = bdx, p_value = stats::pchisq(bdx, K - 1, lower.tail = FALSE))),
                      recommended = list(test = mh$method, rationale = paste0("Binary exposure and outcome compared within ", K,
                                                                               " strata of ", s_nm, ".")),
                      primary_test = add_null_hypothesis(primary, h0),
                      alternative_tests = list(stratum_or = stratum_or,
                                               crude_or = .bt(comparison = cmp, odds_ratio = cr[1], conf.low = cr[2], conf.high = cr[3]),
                                               breslow_day = bdrow("Breslow-Day test", bdx),
                                               breslow_day_tarone = bdrow("Breslow-Day test with Tarone correction", tar)),
                      effect_size = .bt(name = paste0("Mantel-Haenszel OR (", cmp, ", event = ", ev, ")"), estimate = psi,
                                        magnitude = NA_character_),
                      plot = plt, call = match.call(), subclass = "stratified")
  out$interpretation <- make_report(out, alpha)
  out
}
