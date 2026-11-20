cohens_d_independent <- function(data, outcome, group) {
  groups <- assert_two_groups(data, group)
  x <- data[[outcome]][data[[group]] == groups[1]]
  y <- data[[outcome]][data[[group]] == groups[2]]
  nx <- sum(!is.na(x)); ny <- sum(!is.na(y))
  pooled <- sqrt(((nx - 1) * stats::var(x, na.rm = TRUE) + (ny - 1) * stats::var(y, na.rm = TRUE)) / (nx + ny - 2))
  d <- (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) / pooled
  # Hedges g: exact small-sample factor J(df); CIs invert the noncentral t of the pooled t statistic
  # (as effectsize::cohens_d()/hedges_g()).
  df <- nx + ny - 2; h <- sqrt(1 / nx + 1 / ny)
  j <- exp(lgamma(df / 2) - log(sqrt(df / 2)) - lgamma((df - 1) / 2))
  ci <- ncp_ci(d / h, df) * h
  .bt(name = c("Cohen's d", "Hedges' g"), estimate = c(d, d * j), magnitude = magnitude_cohens_d(c(d, d * j)),
      conf.low = c(ci[1], ci[1] * j), conf.high = c(ci[2], ci[2] * j),
      comparison = paste(groups[1], "vs", groups[2]))
}

# CI for the noncentrality parameter of a t statistic with df degrees of freedom.
ncp_ci <- function(t, df, level = 0.95) {
  a <- (1 - level) / 2
  f <- function(p) stats::uniroot(function(ncp) stats::pt(t, df, ncp) - p, t + c(-5, 5), extendInt = "yes",
                                  tol = 1e-10)$root
  c(f(1 - a), f(a))
}

cohens_d_one_sample <- function(x, mu = 0) {
  d <- (mean(x, na.rm = TRUE) - mu) / stats::sd(x, na.rm = TRUE)
  .bt(name = "Cohen's d", estimate = d, magnitude = magnitude_cohens_d(d))
}

cohens_d_paired <- function(before, after) {
  diff <- after - before
  d <- mean(diff, na.rm = TRUE) / stats::sd(diff, na.rm = TRUE)
  .bt(name = "Cohen's dz", estimate = d, magnitude = magnitude_cohens_d(d))
}

# Matched-pairs rank-biserial r_c = (T+ - T-) / (T+ + T-) on the ranks of the non-zero |differences|.
rank_biserial_paired <- function(before, after) {
  d <- after - before; d <- d[!is.na(d) & d != 0]
  rk <- rank(abs(d)); r <- (sum(rk[d > 0]) - sum(rk[d < 0])) / sum(rk)
  .bt(name = "Matched-pairs rank-biserial r", estimate = r, magnitude = magnitude_cramers_v(abs(r)))
}

eta_squared_aov <- function(aov_fit, tab = NULL) {
  # `tab` lets callers pass a pre-computed Type II/III car::Anova() table so
  # eta squared uses the same sums of squares as the reported test, instead
  # of always recomputing Type I SS via stats::anova() regardless of the
  # requested type. Type III tables carry a leading "(Intercept)" row, which
  # is not a treatment effect and must be excluded from both the numerator
  # and the eta-squared denominator.
  if (is.null(tab)) tab <- stats::anova(aov_fit)
  tab <- tab[rownames(tab) != "(Intercept)", , drop = FALSE]
  ss <- tab[["Sum Sq"]]
  eta <- ss[seq_len(length(ss) - 1)] / sum(ss, na.rm = TRUE)
  .bt(term = rownames(tab)[seq_along(eta)], name = "eta squared", estimate = eta, magnitude = magnitude_eta2(eta))
}

cramers_v <- function(tab, biascorrect = FALSE) {
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  if (isTRUE(biascorrect)) {
    # Bergsma (2013): shrink phi^2 by its expected value under independence
    # and use bias-corrected row and column counts.
    r <- nrow(tab); k <- ncol(tab)
    phi2c <- max(0, unname(chi$statistic) / n - (k - 1) * (r - 1) / (n - 1))
    rc <- r - (r - 1)^2 / (n - 1)
    kc <- k - (k - 1)^2 / (n - 1)
    v <- sqrt(phi2c / min(kc - 1, rc - 1))
    return(.bt(name = "Cramer's V (bias-corrected)", estimate = v, magnitude = magnitude_cramers_v(v)))
  }
  k <- min(nrow(tab), ncol(tab))
  v <- sqrt(unname(chi$statistic) / (n * (k - 1)))
  .bt(name = "Cramer's V", estimate = v, magnitude = magnitude_cramers_v(v))
}

rank_biserial_two_groups <- function(data, outcome, group, level = 0.95) {
  groups <- assert_two_groups(data, group)
  x <- data[[outcome]][data[[group]] == groups[1]]; x <- x[!is.na(x)]
  y <- data[[outcome]][data[[group]] == groups[2]]; y <- y[!is.na(y)]
  # theta = P(X > Y) + P(X = Y)/2 = W/(n1 n2) for the first group; r = 2 theta - 1 is positive when the
  # first group tends to have larger values (same direction as Cohen's d). The CI is DeLong's
  # (structural components) Wald interval for theta, mapped to r and clipped to [0, 1].
  psi <- outer(x, y, ">") + 0.5 * outer(x, y, "==")
  th <- mean(psi)
  se <- sqrt(stats::var(rowMeans(psi)) / length(x) + stats::var(colMeans(psi)) / length(y))
  ci <- pmin(1, pmax(0, th + c(-1, 1) * stats::qnorm(1 - (1 - level) / 2) * se))
  r <- 2 * th - 1
  ps <- paste0("Probability of superiority P(", groups[1], " > ", groups[2], ")")
  .bt(name = c("Rank-biserial correlation", ps),
      estimate = c(r, th), magnitude = c(magnitude_cramers_v(abs(r)), NA_character_),
      conf.low = c(2 * ci[1] - 1, ci[1]), conf.high = c(2 * ci[2] - 1, ci[2]),
      comparison = paste(groups[1], "vs", groups[2]))
}

r_squared_lm <- function(fit) {
  g <- broom::glance(fit)
  .bt(name = "R squared", estimate = g$r.squared, magnitude = magnitude_eta2(g$r.squared))
}

adjusted_r_squared_lm <- function(fit) {
  g <- broom::glance(fit)
  .bt(name = "Adjusted R squared", estimate = g$adj.r.squared, magnitude = magnitude_eta2(g$adj.r.squared))
}

mcfadden_r2_glm <- function(fit, null_fit) {
  r2 <- 1 - as.numeric(stats::logLik(fit)) / as.numeric(stats::logLik(null_fit))
  .bt(name = "McFadden's pseudo R squared", estimate = r2, magnitude = magnitude_eta2(r2))
}
