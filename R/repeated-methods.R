repeated_anova_test <- function(long, outcome, within, id) {
  formula <- stats::as.formula(paste(outcome, "~", within, "+ Error(", id, "/", within, ")"))
  fit <- stats::aov(formula, data = long)
  tab <- purrr::map(fit, function(stratum) {
    tryCatch(summary(stratum), error = function(e) NULL)
  })
  tab <- purrr::compact(tab)
  candidates <- purrr::map(tab, 1)
  has_test <- purrr::map_lgl(candidates, function(x) {
    "Pr(>F)" %in% names(x) && any(!is.na(x[["Pr(>F)"]]))
  })
  within_tab <- candidates[[which(has_test)[1]]]
  row <- which(!is.na(within_tab[["Pr(>F)"]]))[1]
  p <- within_tab[["Pr(>F)"]][row]
  statistic <- within_tab[["F value"]][row]
  df1 <- within_tab[["Df"]][row]
  df2 <- within_tab[["Df"]][row + 1]
  tibble::tibble(
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
  tibble::tibble(name = "eta squared", estimate = eta, magnitude = magnitude_eta2(eta))
}

friedman_effect_size <- function(friedman, n, k) {
  w <- unname(friedman$statistic) / (n * (k - 1))
  tibble::tibble(name = "Kendall's W", estimate = w, magnitude = magnitude_cramers_v(w))
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
  tibble::tibble(name = "Cochran Q Kendall's W", estimate = w, magnitude = magnitude_cramers_v(w))
}
