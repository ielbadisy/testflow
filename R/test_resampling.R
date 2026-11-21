RS_STATS <- list(mean = mean, median = stats::median, sd = stats::sd, var = stats::var,
                 trimmed = function(x) mean(x, trim = 0.1))

# Outcome and optional group from `y ~ 1` or `y ~ group`.
rs_vars <- function(formula) {
  f <- stats::as.formula(paste(deparse(formula), collapse = " "))
  list(y = all.vars(f[[2]])[1], g = all.vars(f[[3]])[1], cols = all.vars(f)[1:2][!is.na(all.vars(f)[1:2])])
}

rs_label <- function(statistic) if (is.function(statistic)) "statistic" else statistic

# Quantiles of bootstrap replicates, interpolated on the normal scale as in boot:::norm.inter().
rs_quantile <- function(t, q) {
  t <- sort(t); R <- length(t); rq <- (R + 1) * q; k <- pmin(pmax(floor(rq), 1), R - 1)
  f <- function(j) stats::qnorm(j / (R + 1))
  ifelse(abs(rq - round(rq)) < 1e-10 & round(rq) >= 1 & round(rq) <= R, t[pmin(pmax(round(rq), 1), R)],
         t[k] + (stats::qnorm(q) - f(k)) / (f(k + 1) - f(k)) * (t[k + 1] - t[k]))
}

#' Bootstrap confidence intervals for one- and two-sample statistics
#'
#' Nonparametric bootstrap of a statistic of one sample (`y ~ 1`) or of the
#' difference (or ratio) of a statistic between two groups (`y ~ group`),
#' with normal, percentile, basic and BCa intervals.
#'
#' @param formula `y ~ 1` for one sample or `y ~ group` (two levels) for two
#'   samples.
#' @param data A data frame.
#' @param statistic `"mean"`, `"median"`, `"sd"`, `"var"`, `"trimmed"` (10%
#'   trimmed mean), or a function of a numeric vector returning a scalar.
#'   With two groups the estimand is `statistic(second level) -
#'   statistic(first level)` (or their ratio, see `contrast`).
#' @param type Interval reported as the primary result: `"bca"` (default),
#'   `"percentile"`, `"basic"` or `"normal"`. All four are always in
#'   `$alternative_tests$intervals`.
#' @param B Number of bootstrap resamples.
#' @param contrast `"difference"` (default) or `"ratio"`, for two groups.
#' @param alpha One minus the confidence level.
#' @param seed Optional seed for reproducible resampling.
#' @param plot Logical; include a histogram of the replicates.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_bootstrap`. `$primary_test`
#'   holds the estimate, bootstrap SE (`std.error`), `bias`, the chosen
#'   interval (`conf.low`, `conf.high`, `conf.level`) and `B`;
#'   `$alternative_tests$intervals` has one row per interval type and
#'   `$alternative_tests$bca` the bias correction `z0` and acceleration `a`.
#'   The replicates are in `$replicates`. There is no p-value.
#' @details
#' Two groups are resampled separately, so the group sizes stay fixed (a
#' stratified bootstrap). With replicates \eqn{\hat\theta^*_b}, the basic
#' interval is \eqn{(2\hat\theta - \hat\theta^*_{(1-\alpha/2)},
#' 2\hat\theta - \hat\theta^*_{(\alpha/2)})} and the BCa interval reads the
#' replicates at \eqn{\Phi(z_0 + (z_0 + z_q)/(1 - a(z_0 + z_q)))}, with
#' \eqn{z_0 = \Phi^{-1}\{\#(\hat\theta^* < \hat\theta)/B\}} and the
#' acceleration \eqn{a} from the jackknife (within groups). Quantiles are
#' interpolated on the normal scale. These are the formulas of
#' `boot::boot.ci()` with `L = boot::empinf(type = "jack")`.
#' @references
#' Efron B, Tibshirani RJ. *An Introduction to the Bootstrap*. Chapman &
#' Hall; 1993.
#'
#' Davison AC, Hinkley DV. *Bootstrap Methods and their Application*.
#' Cambridge University Press; 1997.
#' @examples
#' d <- data.frame(y = c(rexp(30), rexp(30, 0.5)), g = rep(c("a", "b"), each = 30))
#' testbootstrap(y ~ g, d, statistic = "median", B = 999, seed = 1, plot = FALSE)
#' @export
testbootstrap <- function(formula, data, statistic = "mean", type = c("bca", "percentile", "basic", "normal"),
                          B = 1999, contrast = c("difference", "ratio"), alpha = 0.05, seed = NULL,
                          plot = TRUE, na.rm = TRUE) {
  type <- match.arg(type); contrast <- match.arg(contrast)
  v <- rs_vars(formula); lab <- rs_label(statistic)
  f <- if (is.function(statistic)) statistic else RS_STATS[[match.arg(statistic, names(RS_STATS))]]
  df <- drop_missing(data, v$cols, na.rm = na.rm)
  y <- df[[v$y]]
  op <- if (contrast == "ratio") `/` else `-`
  if (!is.null(seed)) set.seed(seed)
  if (is.na(v$g)) {
    th <- f(y); n <- length(y)
    t <- vapply(seq_len(B), function(b) f(y[sample.int(n, n, TRUE)]), numeric(1))
    jk <- vapply(seq_len(n), function(i) f(y[-i]), numeric(1))
    L <- (n - 1) * (th - jk)
    name <- paste0(lab, " of ", v$y)
  } else {
    lv <- assert_two_groups(df, v$g); g <- as.character(df[[v$g]])
    y1 <- y[g == lv[1]]; y2 <- y[g == lv[2]]; n1 <- length(y1); n2 <- length(y2)
    th <- op(f(y2), f(y1))
    t <- vapply(seq_len(B), function(b) op(f(y2[sample.int(n2, n2, TRUE)]), f(y1[sample.int(n1, n1, TRUE)])),
                numeric(1))
    L <- c((n1 - 1) * (th - vapply(seq_len(n1), function(i) op(f(y2), f(y1[-i])), numeric(1))),
           (n2 - 1) * (th - vapply(seq_len(n2), function(i) op(f(y2[-i]), f(y1)), numeric(1))))
    name <- paste0(if (contrast == "ratio") "Ratio of " else "Difference in ", lab, " (", lv[2], " vs ", lv[1], ")")
  }
  q <- c(alpha / 2, 1 - alpha / 2); zq <- stats::qnorm(q)
  se <- stats::sd(t); bias <- mean(t) - th
  z0 <- stats::qnorm(mean(t < th)); a <- sum(L^3) / (6 * sum(L^2)^1.5)
  perc <- rs_quantile(t, q)
  ivl <- rbind(normal = th - bias + zq * se, percentile = perc, basic = 2 * th - rev(perc),
               bca = rs_quantile(t, stats::pnorm(z0 + (z0 + zq) / (1 - a * (z0 + zq)))))
  intervals <- .bt(type = rownames(ivl), conf.low = ivl[, 1], conf.high = ivl[, 2])
  label <- paste0("Bootstrap ", c(bca = "BCa", percentile = "percentile", basic = "basic", normal = "normal")[[type]],
                  " interval")
  primary <- .bt(method = label, estimate = th, std.error = se, bias = bias, conf.low = ivl[type, 1],
                 conf.high = ivl[type, 2], conf.level = 1 - alpha, B = B)
  plt <- if (plot) {
    ggplot2::ggplot(data.frame(t = t), ggplot2::aes(x = .data$t)) +
      ggplot2::geom_histogram(bins = 40, fill = "grey75", colour = "white") +
      ggplot2::geom_vline(xintercept = th, colour = TFCOL[["a2"]]) +
      ggplot2::geom_vline(xintercept = ivl[type, ], linetype = "dashed") +
      ggplot2::labs(subtitle = label, x = name, y = "Replicates") +
      tftheme()
  } else NULL
  out <- new_testflow("bootstrap", if (is.na(v$g)) "one-sample bootstrap" else "two-sample bootstrap", v$y,
                      if (is.na(v$g)) NULL else v$g, data = df,
                      descriptives = descriptives_numeric(df, v$y, if (is.na(v$g)) NULL else v$g),
                      assumptions = assumption_checks(check_independence_note(
                        if (is.na(v$g)) "Observations are iid." else "Observations are iid within each group; groups are resampled separately.")),
                      recommended = list(test = label, rationale = "Resampling interval; no distributional form is assumed for the statistic."),
                      primary_test = primary, alternative_tests = list(intervals = intervals, bca = .bt(z0 = z0, a = a)),
                      effect_size = .bt(name = name, estimate = th, magnitude = NA_character_), plot = plt,
                      call = match.call(), subclass = "bootstrap")
  out$replicates <- t
  out$interpretation <- paste0("The ", name, " was ", format_stat(th), " (bootstrap SE ", format_stat(se), ", B = ", B,
                               "). The ", primary_ci_level_text(out), " ", sub("Bootstrap ", "", label), " was [",
                               format_stat(ivl[type, 1]), ", ", format_stat(ivl[type, 2]), "].")
  out
}

# Exact relabellings for k groups: each column of the result is one assignment of group codes.
rs_enumerate <- function(sizes) {
  n <- sum(sizes)
  if (length(sizes) == 1) return(matrix(1L, n, 1))
  rest <- rs_enumerate(sizes[-1])
  idx <- utils::combn(n, sizes[1])
  do.call(cbind, lapply(seq_len(ncol(idx)), function(j) {
    m <- matrix(1L, n, ncol(rest)); m[-idx[, j], ] <- rest + 1L; m
  }))
}

#' Permutation test for two or more independent groups
#'
#' Exact (full enumeration) or Monte Carlo permutation test of a difference
#' between groups, under exchangeability of the outcome across groups.
#'
#' @param formula `y ~ group`.
#' @param data A data frame.
#' @param statistic For two groups: `"mean"` (difference in means, default),
#'   `"median"` (difference in medians) or `"t"` (Welch t statistic), each
#'   as second level minus first level; for more than two groups the
#'   one-way ANOVA F statistic is used. A function `function(y, g)`
#'   returning a scalar is also accepted.
#' @param alternative `"two.sided"`, `"greater"` or `"less"` (two groups only).
#' @param B Number of Monte Carlo relabellings.
#' @param exact `NULL` (default) enumerates all relabellings when there are at
#'   most `1e5`, otherwise uses `B` random ones; `TRUE` or `FALSE` forces
#'   the choice.
#' @param alpha Significance level.
#' @param seed Optional seed for the Monte Carlo relabellings.
#' @param plot Logical; include a histogram of the permutation distribution.
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_permutation`.
#'   `$primary_test` holds the observed `statistic` and the permutation
#'   `p.value`; `$alternative_tests$permutation` gives `exact`, the number of
#'   relabellings and the count at least as extreme. The permutation
#'   distribution is in `$replicates`.
#' @details
#' The exact p-value is the share of all \eqn{n!/\prod n_j!} relabellings
#' whose statistic is at least as extreme as the observed one. The Monte
#' Carlo p-value counts the observed labelling, \eqn{(1 + \#\{|T^*| \ge
#' |t|\})/(B + 1)} (Phipson and Smyth 2010), so it is never 0. The two-sample
#' mean test and the k-sample F test give the same p-values as
#' `coin::oneway_test()`.
#' @references
#' Phipson B, Smyth GK. Permutation p-values should never be zero.
#' *Statistical Applications in Genetics and Molecular Biology*. 2010;9:39.
#'
#' Good PI. *Permutation, Parametric and Bootstrap Tests of Hypotheses*
#' (3rd ed.). Springer; 2005.
#' @examples
#' d <- data.frame(y = c(rnorm(8), rnorm(8, 1)), g = rep(c("a", "b"), each = 8))
#' testpermutation(y ~ g, d, plot = FALSE)
#' @export
testpermutation <- function(formula, data, statistic = "mean", alternative = c("two.sided", "greater", "less"),
                            B = 9999, exact = NULL, alpha = 0.05, seed = NULL, plot = TRUE, na.rm = TRUE) {
  alternative <- match.arg(alternative)
  v <- rs_vars(formula)
  df <- drop_missing(data, v$cols, na.rm = na.rm)
  y <- df[[v$y]]; g <- factor(df[[v$g]]); g <- droplevels(g); lv <- levels(g); k <- length(lv)
  if (k < 2) stop("`", v$g, "` must have at least two groups.", call. = FALSE)
  if (k > 2) alternative <- "greater"
  lab <- if (is.function(statistic)) "statistic" else if (k > 2) "F" else match.arg(statistic, c("mean", "median", "t"))
  fs <- if (is.function(statistic)) statistic else if (k > 2) function(y, g) {
    m <- tapply(y, g, mean); n <- tabulate(g, k)
    sum(n * (m - mean(y))^2) / (k - 1) / (sum((y - m[g])^2) / (length(y) - k))
  } else switch(lab, mean = function(y, g) mean(y[g == 2]) - mean(y[g == 1]),
                median = function(y, g) stats::median(y[g == 2]) - stats::median(y[g == 1]),
                t = function(y, g) unname(stats::t.test(y[g == 2], y[g == 1])$statistic))
  gi <- as.integer(g); sizes <- tabulate(gi, k)
  nperm <- exp(lfactorial(length(y)) - sum(lfactorial(sizes)))
  exact <- exact %||% (nperm <= 1e5)
  tobs <- fs(y, gi)
  if (exact) {
    t <- apply(rs_enumerate(sizes), 2, function(gg) fs(y, gg))
  } else {
    if (!is.null(seed)) set.seed(seed)
    t <- vapply(seq_len(B), function(b) fs(y, sample(gi)), numeric(1))
  }
  tol <- 1e-9 * max(1, abs(tobs))
  hit <- switch(alternative, two.sided = abs(t) >= abs(tobs) - tol, greater = t >= tobs - tol, less = t <= tobs + tol)
  p <- if (exact) mean(hit) else (1 + sum(hit)) / (B + 1)
  slab <- c(mean = "difference in means", median = "difference in medians", t = "Welch t", F = "F",
            statistic = "statistic")[[lab]]
  method <- paste0(if (exact) "Exact" else "Monte Carlo", " permutation test (", slab, ")")
  name <- if (k == 2 && lab %in% c("mean", "median")) paste0("Difference in ", lab, "s (", lv[2], " - ", lv[1], ")") else
    paste0("Observed ", slab)
  primary <- .bt(method = method, estimate = tobs, statistic = tobs, p.value = p, alternative = alternative)
  h0 <- paste0("H0: the distribution of ", v$y, " is the same in every level of ", v$g, " (exchangeable labels).")
  plt <- if (plot) {
    ggplot2::ggplot(data.frame(t = t), ggplot2::aes(x = .data$t)) +
      ggplot2::geom_histogram(bins = 40, fill = "grey75", colour = "white") +
      ggplot2::geom_vline(xintercept = tobs, colour = TFCOL[["a2"]]) +
      ggplot2::labs(subtitle = paste0(method, ", ", pstr(p)), x = slab,
                    y = "Relabellings") +
      tftheme()
  } else NULL
  out <- new_testflow("permutation", if (k == 2) "two independent groups" else "more than two independent groups",
                      v$y, v$g, data = df, descriptives = descriptives_numeric(df, v$y, v$g),
                      assumptions = assumption_checks(assumption_check("Exchangeability", "assumed",
                        "Under H0 group labels are exchangeable (randomisation); unequal spreads can distort a raw mean-difference test.")),
                      recommended = list(test = method, rationale = if (exact) "All relabellings were enumerated." else
                        paste0("Too many relabellings to enumerate; ", B, " random relabellings were drawn.")),
                      primary_test = add_null_hypothesis(primary, h0),
                      alternative_tests = list(permutation = .bt(exact = exact, n_permutations = if (exact) length(t) else B,
                                                                 n_extreme = sum(hit), total_relabellings = nperm)),
                      effect_size = .bt(name = name, estimate = tobs, magnitude = NA_character_), plot = plt,
                      call = match.call(), subclass = "permutation")
  out$replicates <- t
  out$interpretation <- make_report(out, alpha)
  out
}
