POSTHOC_AUTO <- c("One-way ANOVA" = "tukey", "Welch ANOVA" = "welch", "Kruskal-Wallis test" = "wilcoxon")

posthoc_groups <- function(data, outcome, group, recommendation, alpha = 0.05, method = "auto",
                           p.adjust.method = "BH") {
  if (method == "auto") method <- POSTHOC_AUTO[[recommendation]]
  y <- data[[outcome]]; g <- data[[group]]
  adj <- paste0(" (", p.adjust.method, "-adjusted)")
  switch(method,
    tukey = {
      res <- stats::TukeyHSD(stats::aov(tf_formula(outcome, group), data = data))
      # TukeyHSD rows follow lower.tri() order, the same as combn(): "level j - level i".
      tb <- res[[1]]; pp <- ph_pairs(g)
      ph_out("Tukey HSD", res, "Tukey", ph_table(pp$lv[pp$i], pp$lv[pp$j], tb[, "diff"], NA, tb[, "lwr"],
                                                  tb[, "upr"], NA, tb[, "p adj"], "Tukey"))
    },
    gameshowell = {
      tb <- games_howell(y, g)
      ph_out("Games-Howell", ph_htest(tb, "Games-Howell test", "Games-Howell"), "Games-Howell", tb)
    },
    dunn = {
      tb <- dunn_table(y, g, p.adjust.method)
      ph_out(paste0("Dunn test", adj), ph_htest(tb, "Dunn test (pooled ranks)", p.adjust.method),
             p.adjust.method, tb)
    },
    conover = {
      tb <- conover_table(y, g, p.adjust.method)
      ph_out(paste0("Conover-Iman test", adj), ph_htest(tb, "Conover-Iman test (pooled ranks)", p.adjust.method),
             p.adjust.method, tb)
    },
    dscf = {
      tb <- dscf_table(y, g)
      ph_out("Dwass-Steel-Critchlow-Fligner test", ph_htest(tb, "Dwass-Steel-Critchlow-Fligner test", "DSCF"),
             "DSCF", tb)
    },
    welch = {
      if (any(table(g) < 2)) stop("Pairwise Welch t tests need at least 2 observations per group.", call. = FALSE)
      ph_pairwise(stats::pairwise.t.test(y, g, pool.sd = FALSE, p.adjust.method = p.adjust.method),
                  paste0("Pairwise Welch t-tests", adj))
    },
    wilcoxon = ph_pairwise(stats::pairwise.wilcox.test(y, g, p.adjust.method = p.adjust.method, exact = FALSE),
                           paste0("Pairwise Wilcoxon rank-sum tests", adj))
  )
}

ph_out <- function(method, result, padj, table) list(method = method, result = result, p.adjust.method = padj,
                                                      table = table)

ph_table <- function(g1, g2, est, stat, lo, hi, p, padj, m) {
  .bt(group1 = unname(g1), group2 = unname(g2), estimate = unname(est), statistic = unname(stat),
      conf.low = unname(lo), conf.high = unname(hi), p = unname(p), p.adj = unname(padj), p.adjust.method = m)
}

ph_pairwise <- function(pw, method) {
  pm <- pw$p.value; ix <- which(!is.na(pm), arr.ind = TRUE)
  ph_out(method, pw, pw$p.adjust.method, ph_table(colnames(pm)[ix[, 2]], rownames(pm)[ix[, 1]], NA, NA, NA, NA,
                                                  NA, pm[ix], pw$p.adjust.method))
}

# pairwise.htest wrapper (lower-triangular p matrix) so `result` has the same shape for every method.
ph_htest <- function(tb, method, padj) {
  lv <- unique(c(tb$group1, tb$group2))
  pm <- matrix(NA_real_, length(lv) - 1, length(lv) - 1, dimnames = list(lv[-1], lv[-length(lv)]))
  pm[cbind(tb$group2, tb$group1)] <- tb$p.adj
  structure(list(method = method, data.name = "data[[outcome]] and data[[group]]", p.value = pm,
                 p.adjust.method = padj), class = "pairwise.htest")
}

ph_pairs <- function(g) {
  lv <- levels(droplevels(as.factor(g))); cb <- utils::combn(length(lv), 2)
  list(lv = lv, i = cb[1, ], j = cb[2, ])
}

# Games-Howell: Welch t and df per pair, studentized-range p and 95% CI (estimate = mean2 - mean1).
games_howell <- function(y, g, level = 0.95) {
  pp <- ph_pairs(g); g <- factor(g, pp$lv)
  m <- tapply(y, g, mean); v <- tapply(y, g, stats::var); n <- tapply(y, g, length); k <- length(pp$lv)
  i <- pp$i; j <- pp$j; w <- v / n
  est <- m[j] - m[i]; se <- sqrt(w[i] + w[j])
  df <- (w[i] + w[j])^2 / (w[i]^2 / (n[i] - 1) + w[j]^2 / (n[j] - 1))
  t <- est / se; q <- stats::qtukey(level, k, df) / sqrt(2)
  p <- stats::ptukey(abs(t) * sqrt(2), k, df, lower.tail = FALSE)
  ph_table(pp$lv[i], pp$lv[j], est, t, est - q * se, est + q * se, p, p, "Games-Howell")
}

# Dunn (1964) z on ranks pooled over all groups, with tie correction; estimate = mean rank 2 - mean rank 1.
dunn_table <- function(y, g, p.adjust.method = "BH") {
  pp <- ph_pairs(g); g <- factor(g, pp$lv)
  r <- rank(y); N <- length(y); tt <- table(y)
  s2 <- N * (N + 1) / 12 - sum(tt^3 - tt) / (12 * (N - 1))
  rb <- tapply(r, g, mean); n <- tapply(r, g, length); i <- pp$i; j <- pp$j
  est <- rb[j] - rb[i]; z <- est / sqrt(s2 * (1 / n[i] + 1 / n[j])); p <- 2 * stats::pnorm(-abs(z))
  ph_table(pp$lv[i], pp$lv[j], est, z, NA, NA, p, stats::p.adjust(p, p.adjust.method), p.adjust.method)
}

# Conover-Iman (1979) t on pooled ranks: the rank variance S2 and the tie-corrected H shrink the error term.
conover_table <- function(y, g, p.adjust.method = "BH") {
  pp <- ph_pairs(g); g <- factor(g, pp$lv)
  r <- rank(y); N <- length(y); k <- length(pp$lv)
  rb <- tapply(r, g, mean); n <- tapply(r, g, length); i <- pp$i; j <- pp$j
  s2 <- (sum(r^2) - N * (N + 1)^2 / 4) / (N - 1); h <- (sum(n * rb^2) - N * (N + 1)^2 / 4) / s2
  est <- rb[j] - rb[i]; t <- est / sqrt(s2 * (N - 1 - h) / (N - k) * (1 / n[i] + 1 / n[j]))
  p <- 2 * stats::pt(-abs(t), N - k)
  ph_table(pp$lv[i], pp$lv[j], est, t, NA, NA, p, stats::p.adjust(p, p.adjust.method), p.adjust.method)
}

# DSCF: each pair re-ranked, Mann-Whitney U of group2 standardised with the tie-corrected variance, studentized
# range p with k means and infinite df; estimate = Hodges-Lehmann shift (median of group2 - group1 differences).
dscf_table <- function(y, g) {
  pp <- ph_pairs(g); g <- factor(g, pp$lv); k <- length(pp$lv); ys <- split(y, g)
  one <- function(a, b) {
    m <- length(a); n <- length(b); r <- rank(c(a, b)); s <- m + n; tt <- table(r)
    u <- sum(r[m + seq_len(n)]) - n * (n + 1) / 2
    v <- m * n / (s * (s - 1)) * ((s^3 - s) / 12 - sum(tt^3 - tt) / 12)
    c(stats::median(outer(b, a, "-")), sqrt(2) * (u - m * n / 2) / sqrt(v))
  }
  w <- mapply(function(i, j) one(ys[[i]], ys[[j]]), pp$i, pp$j)
  p <- stats::ptukey(abs(w[2, ]), k, Inf, lower.tail = FALSE)
  ph_table(pp$lv[pp$i], pp$lv[pp$j], w[1, ], w[2, ], NA, NA, p, p, "DSCF")
}

# Friedman-type all-pairs tests on within-subject ranks (rows = subjects, columns = conditions, complete cases).
friedman_pairs <- function(wide, method = c("nemenyi", "conover"), p.adjust.method = "BH") {
  method <- match.arg(method)
  r <- t(apply(as.matrix(wide), 1, rank)); n <- nrow(r); k <- ncol(r); lv <- colnames(r)
  cb <- utils::combn(k, 2); i <- cb[1, ]; j <- cb[2, ]; rb <- colMeans(r); est <- rb[j] - rb[i]
  if (method == "nemenyi") {
    q <- sqrt(2) * est / sqrt(k * (k + 1) / (6 * n)); p <- stats::ptukey(abs(q), k, Inf, lower.tail = FALSE)
    return(.bt(group1 = lv[i], group2 = lv[j], method = "Nemenyi test", estimate = unname(est),
               statistic = unname(q), p = unname(p), p.adj = unname(p), p.adjust.method = "Nemenyi"))
  }
  rs <- colSums(r); t <- (rs[j] - rs[i]) / sqrt(2 * (n * sum(r^2) - sum(rs^2)) / ((n - 1) * (k - 1)))
  p <- 2 * stats::pt(-abs(t), (n - 1) * (k - 1))
  .bt(group1 = lv[i], group2 = lv[j], method = "Conover test", estimate = unname(est), statistic = unname(t),
      p = unname(p), p.adj = stats::p.adjust(unname(p), p.adjust.method), p.adjust.method = p.adjust.method)
}

#' Dunn's test for pairwise comparisons after a Kruskal-Wallis test
#'
#' Pairwise z tests on mean ranks computed from the ranks pooled over all
#' groups (Dunn 1964), with the tie correction of the rank variance. Unlike
#' pairwise Wilcoxon rank-sum tests, each pair is not re-ranked, so the tests
#' use the same ranks as the Kruskal-Wallis test.
#' @param formula A formula such as `outcome ~ group`, or a data frame when using data-first style.
#' @param data A data frame, or an outcome column when using data-first style.
#' @param group Grouping column. Optional when using formula style.
#' @param p.adjust.method Multiplicity adjustment passed to [stats::p.adjust()].
#' @param na.rm Logical; remove missing values.
#' @return A table with one row per pair: `group1`, `group2`, `estimate`
#'   (mean rank of `group2` minus mean rank of `group1`), `statistic` (z),
#'   `conf.low`, `conf.high` (both `NA`), `p` (two-sided), `p.adj` and
#'   `p.adjust.method`.
#' @references
#' Dunn, O. J. (1964). Multiple comparisons using rank sums.
#' \emph{Technometrics}, 6(3), 241-252.
#' @examples
#' testdunn(len ~ supp, data = transform(ToothGrowth, supp = interaction(supp, dose)))
#' @export
testdunn <- function(formula, data, group = NULL, p.adjust.method = "BH", na.rm = TRUE) {
  vars <- resolve_formula_pair(tf_formula_arg(substitute(formula), parent.frame()),
                               tf_formula_arg(substitute(data), parent.frame()),
                               substitute(group), missing(group))
  df <- drop_missing(resolve_data_first_or_formula(formula, data), c(vars$outcome, vars$group), na.rm = na.rm)
  dunn_table(df[[vars$outcome]], df[[vars$group]], p.adjust.method)
}

paired_posthoc_numeric <- function(wide, method = c("t", "wilcox"), p.adjust.method = "BH") {
  method <- match.arg(method)
  measures <- names(wide)
  pairs <- utils::combn(measures, 2, simplify = FALSE)
  out <- .map_dfr(pairs, function(pair) {
    x <- wide[[pair[1]]]
    y <- wide[[pair[2]]]
    test <- if (method == "t") {
      stats::t.test(y, x, paired = TRUE)
    } else {
      stats::wilcox.test(y, x, paired = TRUE, exact = FALSE)
    }
    .bt(
      group1 = pair[1],
      group2 = pair[2],
      method = ifelse(method == "t", "paired t-test", "paired Wilcoxon signed-rank test"),
      statistic = unname(test$statistic[1]),
      p = test$p.value
    )
  })
  out$p.adj <- stats::p.adjust(out$p, method = p.adjust.method)
  out$p.adjust.method <- p.adjust.method
  out
}

pairwise_mcnemar <- function(mat, measure_names, p.adjust.method = "BH") {
  pairs <- utils::combn(seq_along(measure_names), 2, simplify = FALSE)
  out <- .map_dfr(pairs, function(pair) {
    tab <- table(mat[, pair[1]], mat[, pair[2]])
    test <- stats::mcnemar.test(tab)
    .bt(
      group1 = measure_names[pair[1]],
      group2 = measure_names[pair[2]],
      method = "McNemar test",
      statistic = unname(test$statistic[1]),
      p = test$p.value
    )
  })
  out$p.adj <- stats::p.adjust(out$p, method = p.adjust.method)
  out$p.adjust.method <- p.adjust.method
  out
}
