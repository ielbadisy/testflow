#' Equivalence and non-inferiority test (TOST) for two independent groups
#'
#' Two one-sided tests (TOST) of whether the difference between two groups
#' lies inside a prespecified margin. Works for means (t-based) and for
#' proportions (Wald z). Setting one bound to `-Inf` or `Inf` gives a
#' one-sided non-inferiority (or non-superiority) test.
#'
#' @param formula A formula `outcome ~ group` with a two-level group.
#' @param data A data frame.
#' @param lower,upper Equivalence bounds on the difference scale
#'   (mean difference, or risk difference for proportions). The difference
#'   is first group level minus second group level, as in [stats::t.test()].
#'   For non-inferiority where higher is better, give the margin as
#'   `lower = -margin, upper = Inf`.
#' @param type `"auto"` (default) uses `"proportions"` when `success` is
#'   given or the outcome is logical, and `"means"` otherwise.
#' @param success Level of the outcome counted as an event for
#'   `type = "proportions"`.
#' @param var.equal Logical; for means, `FALSE` (default) uses Welch
#'   standard errors and df, `TRUE` the pooled variance.
#' @param alpha Significance level of each one-sided test. The reported CI
#'   has level `1 - 2 * alpha` for equivalence (90% at 0.05) and
#'   `1 - alpha` one-sided for non-inferiority.
#' @param plot Logical; include a ggplot object (CI against the margins).
#' @param na.rm Logical; remove missing values.
#' @return A `testflow` object with class `testflow_equivalence`. The
#' primary result is the binding one-sided test (larger p-value), with the
#' difference CI. `$alternative_tests$one_sided` holds both one-sided tests.
#' @details
#' With estimated difference \eqn{\hat\Delta} and standard error \eqn{SE},
#' the two one-sided tests are
#' \deqn{t_L = (\hat\Delta - \Delta_L)/SE, \quad H_{0L}: \Delta \le \Delta_L}
#' \deqn{t_U = (\hat\Delta - \Delta_U)/SE, \quad H_{0U}: \Delta \ge \Delta_U}
#' and equivalence is declared when both are rejected at level
#' \eqn{\alpha}, so \eqn{p_{TOST} = \max(p_L, p_U)}. This is the same as the
#' \eqn{1 - 2\alpha} CI lying inside \eqn{(\Delta_L, \Delta_U)}. For means,
#' \eqn{SE} and df are Welch's (or pooled with `var.equal = TRUE`). For
#' proportions, \eqn{SE = \sqrt{p_1(1-p_1)/n_1 + p_2(1-p_2)/n_2}} with
#' normal quantiles (Wald); it can be poor for small samples or proportions
#' near 0 or 1.
#' @references
#' Schuirmann DJ. A comparison of the two one-sided tests procedure and the
#' power approach for assessing the equivalence of average bioavailability.
#' *Journal of Pharmacokinetics and Biopharmaceutics*. 1987;15(6):657-680.
#'
#' Walker E, Nowacki AS. Understanding equivalence and noninferiority
#' testing. *Journal of General Internal Medicine*. 2011;26(2):192-196.
#' @examples
#' set.seed(2)
#' d <- data.frame(sbp = c(rnorm(60, 130, 12), rnorm(60, 131, 12)),
#'                 arm = rep(c("new", "standard"), each = 60))
#' testequivalence(sbp ~ arm, d, lower = -5, upper = 5, plot = FALSE)
#' # non-inferiority: new arm no worse than 5 mmHg higher
#' testequivalence(sbp ~ arm, d, lower = -Inf, upper = 5, plot = FALSE)
#' @export
testequivalence <- function(formula, data, lower, upper, type = c("auto", "means", "proportions"), success = NULL, var.equal = FALSE, alpha = 0.05, plot = TRUE, na.rm = TRUE) {
  type <- match.arg(type)
  if (!inherits(formula, "formula")) stop("`formula` must be a formula such as outcome ~ group.", call. = FALSE)
  if (missing(lower) || missing(upper)) stop("Give both `lower` and `upper` (use -Inf or Inf for non-inferiority).", call. = FALSE)
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != 1 || length(upper) != 1 || lower >= upper) {
    stop("`lower` must be a number smaller than `upper`.", call. = FALSE)
  }
  if (is.infinite(lower) && is.infinite(upper)) stop("At most one bound can be infinite.", call. = FALSE)
  outcome_nm <- all.vars(formula[[2]])[1]
  group_nm <- all.vars(formula[[3]])[1]
  df <- drop_missing(data, c(outcome_nm, group_nm), na.rm = na.rm)
  df[[group_nm]] <- as.factor(df[[group_nm]])
  assert_two_groups(df, group_nm)
  lv <- levels(droplevels(df[[group_nm]]))
  y <- df[[outcome_nm]]
  if (type == "auto") type <- if (!is.null(success) || is.logical(y)) "proportions" else "means"
  g1 <- df[[group_nm]] == lv[1]
  g2 <- df[[group_nm]] == lv[2]
  n1 <- sum(g1)
  n2 <- sum(g2)

  if (type == "means") {
    if (!is.numeric(y)) stop("A numeric outcome is needed for type = \"means\".", call. = FALSE)
    m1 <- mean(y[g1]); m2 <- mean(y[g2])
    v1 <- stats::var(y[g1]); v2 <- stats::var(y[g2])
    est <- m1 - m2
    if (var.equal) {
      sp2 <- ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
      se <- sqrt(sp2 * (1 / n1 + 1 / n2))
      dfree <- n1 + n2 - 2
    } else {
      se <- sqrt(v1 / n1 + v2 / n2)
      dfree <- (v1 / n1 + v2 / n2)^2 / ((v1 / n1)^2 / (n1 - 1) + (v2 / n2)^2 / (n2 - 1))
    }
    upper_q <- function(p) stats::qt(p, dfree)
    p_lower_fn <- function(s) stats::pt(s, dfree, lower.tail = FALSE)
    p_upper_fn <- function(s) stats::pt(s, dfree)
    method <- paste0("TOST, ", if (var.equal) "Student" else "Welch", " t (difference in means)")
    est_name <- "Mean difference"
  } else {
    if (is.null(success)) {
      if (!is.logical(y)) stop("Give `success` for type = \"proportions\".", call. = FALSE)
      success <- TRUE
    }
    ev <- y == success
    p1 <- mean(ev[g1]); p2 <- mean(ev[g2])
    est <- p1 - p2
    se <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    if (se == 0) stop("Standard error is zero (all or no events in both groups); the Wald TOST is undefined.", call. = FALSE)
    dfree <- NA_real_
    upper_q <- function(p) stats::qnorm(p)
    p_lower_fn <- function(s) stats::pnorm(s, lower.tail = FALSE)
    p_upper_fn <- function(s) stats::pnorm(s)
    method <- "TOST, Wald z (difference in proportions)"
    est_name <- "Risk difference"
  }

  t_low <- (est - lower) / se
  t_up <- (est - upper) / se
  p_low <- if (is.finite(lower)) p_lower_fn(t_low) else 0
  p_up <- if (is.finite(upper)) p_upper_fn(t_up) else 0
  objective <- if (is.finite(lower) && is.finite(upper)) "equivalence" else "noninferiority"
  q <- upper_q(1 - alpha)
  ci <- c(if (is.finite(lower)) est - q * se else -Inf, if (is.finite(upper)) est + q * se else Inf)
  if (objective == "equivalence") ci <- est + c(-1, 1) * q * se
  binding_low <- p_low >= p_up
  p_tost <- max(p_low, p_up)
  conclusion <- p_tost < alpha

  diff_lab <- paste0(lv[1], " - ", lv[2])
  bounds_txt <- paste0("[", format_stat(lower), ", ", format_stat(upper), "]")
  h0 <- if (objective == "equivalence") {
    paste0("H0: the difference (", diff_lab, ") lies outside ", bounds_txt, " (not equivalent).")
  } else if (is.finite(lower)) {
    paste0("H0: the difference (", diff_lab, ") is at most ", format_stat(lower), " (inferior).")
  } else {
    paste0("H0: the difference (", diff_lab, ") is at least ", format_stat(upper), " (inferior).")
  }
  one_sided <- .bt(
    test = c("lower bound", "upper bound"),
    bound = c(lower, upper),
    statistic = c(if (is.finite(lower)) t_low else NA_real_, if (is.finite(upper)) t_up else NA_real_),
    parameter = dfree,
    p.value = c(if (is.finite(lower)) p_low else NA_real_, if (is.finite(upper)) p_up else NA_real_),
    null_hypothesis = c(paste0("H0: difference <= ", format_stat(lower)), paste0("H0: difference >= ", format_stat(upper)))
  )
  primary <- .bt(
    method = method,
    estimate = est,
    statistic = if (binding_low) t_low else t_up,
    parameter = dfree,
    p.value = p_tost,
    conf.low = ci[1],
    conf.high = ci[2],
    conf.level = if (objective == "equivalence") 1 - 2 * alpha else 1 - alpha,
    null_hypothesis = h0
  )
  effect <- .bt(name = paste0(est_name, " (", diff_lab, ")"), estimate = est, magnitude = NA_character_)
  ci_level <- if (objective == "equivalence") 1 - 2 * alpha else 1 - alpha

  plt <- if (plot) {
    # Thick bar: the interval the decision uses (90% for TOST at alpha = 0.05); thin bar: the usual two-sided CI.
    w <- function(q) c(est - q * se, est + q * se); i1 <- w(q); i2 <- w(upper_q(1 - alpha / 2))
    fb <- c(lower, upper)[is.finite(c(lower, upper))]
    rng <- range(c(i2, fb, 0)); pad <- diff(rng) * 0.08
    ggplot2::ggplot() +
      ggplot2::annotate("rect", xmin = if (is.finite(lower)) lower else rng[1] - pad,
                        xmax = if (is.finite(upper)) upper else rng[2] + pad, ymin = -Inf, ymax = Inf,
                        fill = TFCOL[["a1"]], alpha = 0.08) +
      ggplot2::geom_vline(xintercept = fb, linetype = "dashed", colour = TFCOL[["a1"]]) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
      ggplot2::annotate("segment", x = i2[1], xend = i2[2], y = 1, yend = 1, colour = "grey45", linewidth = 0.5) +
      ggplot2::annotate("segment", x = i1[1], xend = i1[2], y = 1, yend = 1, colour = TFCOL[["ink"]], linewidth = 1.8) +
      ggplot2::annotate("point", x = est, y = 1, size = 3, shape = 21, fill = "white", colour = TFCOL[["ink"]]) +
      ggplot2::scale_y_continuous(limits = c(0.6, 1.4), breaks = NULL) +
      ggplot2::coord_cartesian(xlim = rng + c(-pad, pad)) +
      ggplot2::labs(
        subtitle = paste0(if (objective == "equivalence") "TOST" else "Non-inferiority", ", margin ", bounds_txt,
                          ", ", pstr(p_tost)),
        caption = paste0("Thick: ", format(100 * (1 - 2 * alpha)), "% CI (decision); thin: ",
                         format(100 * (1 - alpha)), "% CI; dashed: margins."),
        x = paste0(est_name, " (", diff_lab, ")"), y = NULL
      ) +
      tftheme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                                 panel.grid.minor = ggplot2::element_blank())
  } else NULL

  out <- new_testflow(
    "equivalence",
    if (objective == "equivalence") "equivalence (TOST)" else "non-inferiority",
    outcome_nm,
    group_nm,
    data = df,
    descriptives = if (type == "means") descriptives_numeric(df, outcome_nm, group_nm) else descriptives_categorical(df, c(outcome_nm, group_nm)),
    assumptions = assumption_checks(
      check_independence_note(),
      assumption_check("Margin fixed in advance", "assumed", "The margin must come from clinical reasoning before seeing the data.")
    ),
    recommended = list(test = method, rationale = paste0("The question is whether the difference lies within ", bounds_txt, ", so absence of a significant difference is not enough; both one-sided nulls are tested.")),
    primary_test = primary,
    alternative_tests = list(one_sided = one_sided),
    effect_size = effect,
    plot = plt,
    call = match.call(),
    subclass = "equivalence"
  )
  out$objective <- objective
  out$margin <- c(lower = lower, upper = upper)
  out$interpretation <- paste0(
    "The ", out$design, " workflow for ", outcome_nm, " estimated a ", tolower(est_name), " (", diff_lab, ") of ",
    format_stat(est), ", ", format(100 * ci_level), "% CI [", format_stat(ci[1]), ", ", format_stat(ci[2]), "], against the margin ",
    bounds_txt, " (", method, ", p = ", format_p(p_tost), "). ",
    if (conclusion) {
      if (objective == "equivalence") "Equivalence was shown at this margin." else "Non-inferiority was shown at this margin."
    } else {
      if (objective == "equivalence") "Equivalence was not shown at this margin." else "Non-inferiority was not shown at this margin."
    },
    " ", h0
  )
  out
}
