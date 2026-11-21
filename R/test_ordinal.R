reg_vars <- function(formula) {
  f <- stats::as.formula(paste(deparse(formula), collapse = " "))
  list(f = f, y = all.vars(f[[2]])[1], x = all.vars(f[[3]]))
}

# Brant (1990) Wald test of parallel slopes from the J - 1 cumulative binary logits.
brant_test <- function(f, df, y) {
  X <- stats::model.matrix(f, df); K <- ncol(X) - 1; yi <- as.integer(y); J <- max(yi)
  fits <- lapply(seq_len(J - 1), function(m) stats::glm.fit(X, as.numeric(yi > m), family = stats::binomial()))
  pis <- vapply(fits, `[[`, numeric(nrow(X)), "fitted.values")
  inv <- lapply(seq_len(J - 1), function(m) solve(crossprod(X, X * (pis[, m] * (1 - pis[, m])))))
  b <- unlist(lapply(fits, function(m) m$coefficients[-1]))
  V <- matrix(0, (J - 1) * K, (J - 1) * K); ix <- function(m) ((m - 1) * K + 1):(m * K)
  for (m in seq_len(J - 1)) for (l in m:(J - 1)) {
    # Diagonal blocks are the glm covariances (IRLS weights of the last iteration), as brant::brant() uses.
    blk <- if (m == l) chol2inv(fits[[m]]$qr$qr[seq_len(K + 1), seq_len(K + 1), drop = FALSE]) else inv[[m]] %*% crossprod(X, X * (pis[, l] - pis[, m] * pis[, l])) %*% inv[[l]]
    V[ix(m), ix(l)] <- blk[-1, -1]; V[ix(l), ix(m)] <- t(blk[-1, -1])
  }
  D <- kronecker(cbind(1, -diag(J - 2)), diag(K))
  wald <- function(D) { d <- D %*% b; c(stat = drop(t(d) %*% solve(D %*% V %*% t(D)) %*% d), df = nrow(D)) }
  rows <- rbind(Omnibus = wald(D), t(vapply(seq_len(K), function(k) wald(D[seq(k, nrow(D), K), , drop = FALSE]), numeric(2))))
  .bt(term = c("Omnibus", colnames(X)[-1]), statistic = rows[, 1], parameter = rows[, 2],
      p.value = stats::pchisq(rows[, 1], rows[, 2], lower.tail = FALSE))
}

# Wald table on the log scale and exponentiated.
wald_tab <- function(term, est, se, alpha, ename, extra = NULL) {
  z <- stats::qnorm(1 - alpha / 2)
  out <- .bt(term = term, estimate = est, std.error = se, statistic = est / se, p.value = 2 * stats::pnorm(-abs(est / se)),
             conf.low = est - z * se, conf.high = est + z * se)
  if (!is.null(extra)) out <- cbind(extra, out)
  ex <- out; ex$estimate <- exp(out$estimate); ex$conf.low <- exp(out$conf.low); ex$conf.high <- exp(out$conf.high)
  names(ex)[names(ex) == "estimate"] <- ename
  list(log = .as_bt(out), exp = .as_bt(ex))
}

#' Ordinal logistic regression (proportional odds) workflow
#'
#' Fits a cumulative-logit proportional odds model with `MASS::polr()`,
#' reports cumulative odds ratios with Wald intervals, the likelihood-ratio
#' test against the null model and a check of the proportional odds
#' assumption.
#'
#' @param formula A formula such as `outcome ~ x1 + x2`. The outcome is used
#'   in its factor level order (lowest first).
#' @param data A data frame.
#' @param potest Proportional odds check: `"brant"` (default; Brant Wald
#'   test, omnibus and per coefficient) or `"lr"` (likelihood-ratio test of
#'   the proportional odds model against the multinomial logit, which
#'   relaxes every slope).
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param units Optional named numeric vector of plotting units for numeric
#'   predictors, for example `c(age = 10)` to plot the odds ratio per 10
#'   years, \eqn{e^{10\beta}}. Only the plot changes; tables stay per unit.
#' @param labels,unitlabels Optional plot labels, see
#'   [testlogisticregression()].
#' @return A `testflow` object with class `testflow_ordinal_regression`.
#'   `$primary_test` is the LR test of all slopes (`statistic`, `parameter`,
#'   `p.value`). `$alternative_tests` has `coefficients` (log odds, Wald
#'   SE, z, p, CI), `odds_ratios` (`odds_ratio`, `conf.low`, `conf.high`),
#'   `thresholds` (cut-points \eqn{\zeta_j}) and `proportional_odds` (Brant
#'   rows `Omnibus` and one per coefficient, or the LR row). `$effect_size`
#'   is McFadden's pseudo \eqn{R^2}; the fit is in `$model`.
#' @details
#' The model is \eqn{\mathrm{logit}\,P(Y \le j) = \zeta_j - x^\top\beta},
#' so \eqn{e^{\beta_k}} is the odds ratio for a higher rather than lower
#' outcome category per unit of \eqn{x_k}, the same at every cut-point. The
#' Brant test compares the slopes of the \eqn{J - 1} binary logits
#' \eqn{P(Y > j)} (Wald, \eqn{(J - 2)K} df). Per-coefficient rows match
#' `brant::brant()`; the omnibus row can differ slightly because brant
#' copies the upper covariance blocks below the diagonal without
#' transposing them, while this function keeps the matrix symmetric.
#' @references
#' McCullagh P. Regression models for ordinal data. *JRSS B*.
#' 1980;42:109-142.
#'
#' Brant R. Assessing proportionality in the proportional odds model for
#' ordinal logistic regression. *Biometrics*. 1990;46:1171-1178.
#' @examples
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   h <- MASS::housing[rep(1:72, MASS::housing$Freq), ]
#'   testordinalregression(Sat ~ Infl + Type, h, plot = FALSE)
#' }
#' @export
testordinalregression <- function(formula, data, potest = c("brant", "lr"), alpha = 0.05, plot = TRUE, na.rm = TRUE,
                                  units = NULL, labels = NULL, unitlabels = NULL) {
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package `MASS` is required.", call. = FALSE)
  potest <- match.arg(potest); v <- reg_vars(formula)
  df <- drop_missing(data, c(v$y, v$x), na.rm = na.rm)
  df[[v$y]] <- droplevels(as.factor(df[[v$y]]))
  if (nlevels(df[[v$y]]) < 3) stop("`", v$y, "` needs at least three ordered levels.", call. = FALSE)
  fit <- MASS::polr(v$f, data = df, Hess = TRUE)
  null_fit <- MASS::polr(stats::update(v$f, . ~ 1), data = df, Hess = TRUE)
  se <- sqrt(diag(stats::vcov(fit))); b <- stats::coef(fit); p <- length(b)
  wt <- wald_tab(names(b), b, se[seq_len(p)], alpha, "odds_ratio")
  thr <- .bt(threshold = names(fit$zeta), estimate = unname(fit$zeta), std.error = unname(se[-seq_len(p)]))
  lr <- 2 * (as.numeric(stats::logLik(fit)) - as.numeric(stats::logLik(null_fit)))
  po <- if (potest == "brant") brant_test(v$f, df, df[[v$y]]) else {
    if (!requireNamespace("nnet", quietly = TRUE)) stop("Package `nnet` is required for potest = \"lr\".", call. = FALSE)
    mn <- nnet::multinom(v$f, data = df, trace = FALSE, maxit = 1000)
    st <- 2 * (as.numeric(stats::logLik(mn)) - as.numeric(stats::logLik(fit))); dd <- mn$edf - fit$edf
    .bt(term = "Omnibus (LR vs multinomial)", statistic = st, parameter = dd, p.value = stats::pchisq(st, dd, lower.tail = FALSE))
  }
  pop <- po$p.value[1]
  po_check <- assumption_check("Proportional odds", if (pop >= alpha) "acceptable" else "not acceptable",
                               if (pop >= alpha) "No evidence that slopes differ across cut-points." else
                                 "Slopes differ across cut-points; consider a partial proportional odds or multinomial model.",
                               method = if (potest == "brant") "Brant test" else "LR test vs multinomial logit",
                               statistic = po$statistic[1], p_value = pop)
  primary <- .bt(method = "Ordinal logistic regression (likelihood-ratio test)", statistic = lr, parameter = p,
                 p.value = stats::pchisq(lr, p, lower.tail = FALSE))
  h0 <- paste0("H0: none of ", paste(v$x, collapse = ", "), " are associated with ", v$y, " (all cumulative log-odds slopes equal 0).")
  plt <- if (plot) ratio_plot(wt$log, paste0("Cumulative odds ratio, ", 100 * (1 - alpha), "% CI (log scale)"), v$x, units,
                              subtitle = paste0("Likelihood-ratio test, ", pstr(primary$p.value)), labels = labels,
                              unitlabels = unitlabels) else NULL
  out <- new_testflow("ordinal_regression", "ordinal logistic regression", v$y, paste(v$x, collapse = ", "), data = df,
                      descriptives = descriptives_categorical(df, v$y),
                      assumptions = assumption_checks(check_independence_note(), po_check, check_multicollinearity(fit)),
                      recommended = list(test = "Ordinal logistic regression (proportional odds)",
                                         rationale = "Ordered categorical outcome; one odds ratio per predictor across all cut-points."),
                      primary_test = add_null_hypothesis(primary, h0),
                      alternative_tests = list(coefficients = wt$log, odds_ratios = wt$exp, thresholds = thr, proportional_odds = po),
                      effect_size = mcfadden_r2_glm(fit, null_fit), plot = plt, call = match.call(), subclass = "ordinal_regression")
  out$model <- fit
  out$interpretation <- make_report(out, alpha)
  out
}

#' Multinomial logistic regression workflow
#'
#' Fits a baseline-category logit model with `nnet::multinom()` and reports
#' relative risk ratios with Wald intervals and the likelihood-ratio test
#' against the null model.
#'
#' @param formula A formula such as `outcome ~ x1 + x2`. The first level of
#'   the outcome is the baseline category.
#' @param data A data frame.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param units Optional named numeric vector of plotting units for numeric
#'   predictors (see [testordinalregression()]).
#' @param labels,unitlabels Optional plot labels, see
#'   [testlogisticregression()].
#' @return A `testflow` object with class `testflow_multinomial_regression`.
#'   `$primary_test` is the LR test of all slopes. `$alternative_tests` has
#'   `coefficients` (log scale, one row per `y.level` and `term`, Wald SE,
#'   z, p, CI) and `relative_risk_ratios` (`rrr`, `conf.low`, `conf.high`).
#'   `$effect_size` is McFadden's pseudo \eqn{R^2}; the fit is in `$model`.
#' @details
#' For each non-baseline category \eqn{j},
#' \eqn{\log\{P(Y = j)/P(Y = 1)\} = \alpha_j + x^\top\beta_j}, and
#' \eqn{e^{\beta_{jk}}} is the relative risk ratio of category \eqn{j}
#' versus the baseline per unit of \eqn{x_k}.
#' @references
#' Agresti A. *Categorical Data Analysis* (3rd ed.). Wiley; 2013. Chapter 8.
#'
#' Venables WN, Ripley BD. *Modern Applied Statistics with S* (4th ed.).
#' Springer; 2002.
#' @examples
#' if (requireNamespace("nnet", quietly = TRUE)) {
#'   testmultinomialregression(Species ~ Sepal.Width, iris, plot = FALSE)
#' }
#' @export
testmultinomialregression <- function(formula, data, alpha = 0.05, plot = TRUE, na.rm = TRUE, units = NULL,
                                      labels = NULL, unitlabels = NULL) {
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Package `nnet` is required.", call. = FALSE)
  v <- reg_vars(formula)
  df <- drop_missing(data, c(v$y, v$x), na.rm = na.rm)
  df[[v$y]] <- droplevels(as.factor(df[[v$y]]))
  if (nlevels(df[[v$y]]) < 3) stop("`", v$y, "` needs at least three levels; use testlogisticregression() for two.", call. = FALSE)
  fit <- nnet::multinom(v$f, data = df, trace = FALSE, Hess = TRUE, maxit = 1000)
  null_fit <- nnet::multinom(stats::update(v$f, . ~ 1), data = df, trace = FALSE, maxit = 1000)
  b <- stats::coef(fit); se <- matrix(sqrt(diag(stats::vcov(fit))), nrow(b), byrow = TRUE)
  wt <- wald_tab(rep(colnames(b), nrow(b)), as.vector(t(b)), as.vector(t(se)), alpha, "rrr",
                 extra = .bt(y.level = rep(rownames(b), each = ncol(b))))
  lr <- 2 * (as.numeric(stats::logLik(fit)) - as.numeric(stats::logLik(null_fit))); dd <- fit$edf - null_fit$edf
  primary <- .bt(method = "Multinomial logistic regression (likelihood-ratio test)", statistic = lr, parameter = dd,
                 p.value = stats::pchisq(lr, dd, lower.tail = FALSE))
  h0 <- paste0("H0: none of ", paste(v$x, collapse = ", "), " are associated with ", v$y, " (all baseline-category slopes equal 0).")
  rrr <- wt$exp[wt$exp$term != "(Intercept)", ]
  plt <- if (plot) ratio_plot(wt$log, paste0("Relative risk ratio vs ", levels(df[[v$y]])[1], ", ", 100 * (1 - alpha),
                                            "% CI (log scale)"), v$x, units,
                              subtitle = paste0("Likelihood-ratio test, ", pstr(primary$p.value)), labels = labels,
                              unitlabels = unitlabels) else NULL
  out <- new_testflow("multinomial_regression", "multinomial logistic regression", v$y, paste(v$x, collapse = ", "), data = df,
                      descriptives = descriptives_categorical(df, v$y),
                      assumptions = assumption_checks(check_independence_note(), assumption_check(
                        "Independence of irrelevant alternatives", "assumed",
                        "Odds between two categories do not depend on the other categories.")),
                      recommended = list(test = "Multinomial logistic regression",
                                         rationale = paste0("Unordered outcome with ", nlevels(df[[v$y]]), " categories; baseline ",
                                                            levels(df[[v$y]])[1], ".")),
                      primary_test = add_null_hypothesis(primary, h0),
                      alternative_tests = list(coefficients = wt$log, relative_risk_ratios = wt$exp),
                      effect_size = mcfadden_r2_glm(fit, null_fit), plot = plt, call = match.call(),
                      subclass = "multinomial_regression")
  out$model <- fit
  out$interpretation <- make_report(out, alpha)
  out
}
