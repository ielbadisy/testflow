cohens_d_independent <- function(data, outcome, group) {
  groups <- assert_two_groups(data, group)
  x <- data[[outcome]][data[[group]] == groups[1]]
  y <- data[[outcome]][data[[group]] == groups[2]]
  nx <- sum(!is.na(x)); ny <- sum(!is.na(y))
  pooled <- sqrt(((nx - 1) * stats::var(x, na.rm = TRUE) + (ny - 1) * stats::var(y, na.rm = TRUE)) / (nx + ny - 2))
  d <- (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) / pooled
  tibble::tibble(name = "Cohen's d", estimate = d, magnitude = magnitude_cohens_d(d))
}

cohens_d_one_sample <- function(x, mu = 0) {
  d <- (mean(x, na.rm = TRUE) - mu) / stats::sd(x, na.rm = TRUE)
  tibble::tibble(name = "Cohen's d", estimate = d, magnitude = magnitude_cohens_d(d))
}

cohens_d_paired <- function(before, after) {
  diff <- after - before
  d <- mean(diff, na.rm = TRUE) / stats::sd(diff, na.rm = TRUE)
  tibble::tibble(name = "Cohen's dz", estimate = d, magnitude = magnitude_cohens_d(d))
}

eta_squared_aov <- function(aov_fit) {
  tab <- stats::anova(aov_fit)
  ss <- tab[["Sum Sq"]]
  eta <- ss[seq_len(length(ss) - 1)] / sum(ss, na.rm = TRUE)
  tibble::tibble(term = rownames(tab)[seq_along(eta)], name = "eta squared", estimate = eta, magnitude = magnitude_eta2(eta))
}

cramers_v <- function(tab) {
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  v <- sqrt(unname(chi$statistic) / (n * (k - 1)))
  tibble::tibble(name = "Cramer's V", estimate = v, magnitude = magnitude_cramers_v(v))
}

rank_biserial_two_groups <- function(data, outcome, group) {
  groups <- assert_two_groups(data, group)
  x <- data[[outcome]][data[[group]] == groups[1]]
  y <- data[[outcome]][data[[group]] == groups[2]]
  wt <- stats::wilcox.test(x, y, exact = FALSE)
  n1 <- sum(!is.na(x)); n2 <- sum(!is.na(y))
  r <- 1 - (2 * unname(wt$statistic)) / (n1 * n2)
  tibble::tibble(name = "Rank-biserial correlation", estimate = r, magnitude = magnitude_cramers_v(abs(r)))
}
