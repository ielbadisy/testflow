repeated_anova_test <- function(long, outcome, within, id) {
  formula <- stats::as.formula(paste(tf_backtick(outcome), "~", tf_backtick(within), "+ Error(", tf_backtick(id), "/", tf_backtick(within), ")"))
  fit <- stats::aov(formula, data = long)
  tab <- basetable::map(fit, function(stratum) {
    tryCatch(summary(stratum), error = function(e) NULL)
  })
  tab <- basetable::compact(tab)
  candidates <- lapply(tab, `[[`, 1)
  has_test <- vapply(candidates, function(x) {
    "Pr(>F)" %in% names(x) && any(!is.na(x[["Pr(>F)"]]))
  }, logical(1))
  within_tab <- candidates[[which(has_test)[1]]]
  row <- which(!is.na(within_tab[["Pr(>F)"]]))[1]
  p <- within_tab[["Pr(>F)"]][row]
  statistic <- within_tab[["F value"]][row]
  df1 <- within_tab[["Df"]][row]
  df2 <- within_tab[["Df"]][row + 1]
  .bt(
    method = "Repeated-measures ANOVA",
    statistic = statistic,
    parameter = df1,
    df.error = df2,
    p.value = p,
    conf.low = NA_real_,
    conf.high = NA_real_
  )
}

repeated_eta_squared <- function(long, outcome, within, id) {
  grand <- mean(long[[outcome]], na.rm = TRUE)
  ss_total <- sum((long[[outcome]] - grand)^2, na.rm = TRUE)
  time_means <- stats::aggregate(long[[outcome]], list(time = long[[within]]), mean, na.rm = TRUE)
  n_by_time <- as.numeric(table(long[[within]]))
  ss_time <- sum(n_by_time * (time_means$x - grand)^2, na.rm = TRUE)
  eta <- ss_time / ss_total
  .bt(name = "eta squared", estimate = eta, magnitude = magnitude_eta2(eta))
}

friedman_effect_size <- function(friedman, n, k) {
  w <- unname(friedman$statistic) / (n * (k - 1))
  .bt(name = "Kendall's W", estimate = w, magnitude = magnitude_cramers_v(w))
}

cochran_q_test <- function(mat) {
  if (!all(mat %in% c(0, 1))) {
    stop("Cochran Q requires binary repeated measures coded as two categories.", call. = FALSE)
  }
  k <- ncol(mat)
  n <- nrow(mat)
  col_totals <- colSums(mat)
  row_totals <- rowSums(mat)
  numerator <- (k - 1) * (k * sum(col_totals^2) - sum(col_totals)^2)
  denominator <- k * sum(row_totals) - sum(row_totals^2)
  q <- if (denominator == 0) NA_real_ else numerator / denominator
  df <- k - 1
  p <- stats::pchisq(q, df = df, lower.tail = FALSE)
  structure(
    list(statistic = c(Q = q), parameter = c(df = df), p.value = p, method = "Cochran Q test", data.name = deparse(substitute(mat))),
    class = "htest"
  )
}

cochran_q_effect_size <- function(test, n, k) {
  w <- unname(test$statistic) / (n * (k - 1))
  .bt(name = "Cochran Q Kendall's W", estimate = w, magnitude = magnitude_cramers_v(w))
}

# Greenhouse-Geisser and Huynh-Feldt epsilon from the complete wide matrix (one row per subject), using
# the residual covariance of lm(M ~ between) on orthonormal within contrasts. HF uses the same
# formula as anova.mlm() and car::Anova(): ((n + 1) p GG - 2) / (p (n - p GG)), n = residual df.
sphericity_eps <- function(M, between = NULL) {
  p <- ncol(M) - 1
  if (p < 2) return(c(gg = 1, hf = 1))
  fit <- if (is.null(between)) stats::lm(M ~ 1) else stats::lm(M ~ factor(between))
  n <- fit$df.residual; C <- stats::contr.poly(ncol(M))
  l <- eigen(crossprod(C, crossprod(stats::residuals(fit)) %*% C), symmetric = TRUE, only.values = TRUE)$values
  gg <- sum(l)^2 / sum(l^2) / p
  c(gg = gg, hf = ((n + 1) * p * gg - 2) / (p * (n - p * gg)))
}

# Adds gg.eps, hf.eps and the corrected p-values (df scaled by epsilon, HF capped at 1) to F-test rows.
add_eps <- function(tb, eps, df1 = tb$parameter, df2 = tb$df.error, rows = TRUE) {
  na <- ifelse(rows, 1, NA)
  tb$gg.eps <- eps[["gg"]] * na; tb$hf.eps <- eps[["hf"]] * na
  pc <- function(e) stats::pf(tb$statistic, df1 * e, df2 * e, lower.tail = FALSE) * na
  tb$p.gg <- pc(eps[["gg"]]); tb$p.hf <- pc(min(1, eps[["hf"]]))
  tb
}
