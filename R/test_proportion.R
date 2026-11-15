#' Test a one-sample proportion
#' @param data A data frame.
#' @param outcome Categorical outcome column.
#' @param success Value counted as success.
#' @param p Reference probability.
#' @param alpha Significance level.
#' @param plot Logical; include a ggplot object.
#' @param na.rm Logical; remove missing values.
#' @param weights,cluster,strata Optional column names (strings) for a
#'   complex survey design: sampling weights, primary sampling units and
#'   strata. When any is given, the primary test becomes design-based (via the
#'   `survey` package): the weighted proportion, its linearisation SE, a
#'   t test of \eqn{H_0: \pi = p} as in `survey::svyttest()`, and a
#'   logit CI from `survey::svyciprop()`. The design effect is reported. The
#'   unweighted exact and approximate tests stay in `alternative_tests`.
#'   Default `NULL` keeps the simple random sample analysis.
#' @return A `testflow` object with class `testflow_proportion`. The object is
#' a list containing the cleaned data, categorical descriptives, assumption
#' checks, recommended exact or approximate one-sample proportion test, primary
#' test result with null hypothesis, alternative exact and approximate results,
#' observed proportion summary, optional `ggplot`, original call, and report
#' text.
#' @references
#' Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or fiducial
#' limits illustrated in the case of the binomial. \emph{Biometrika},
#' 26(4), 404-413.
#'
#' Lumley, T. (2004). Analysis of complex survey samples. \emph{Journal of
#' Statistical Software}, 9(1), 1-19.
#' @export
testproportion <- function(data, outcome, success, p = 0.5, alpha = 0.05, plot = TRUE, na.rm = TRUE,
                           weights = NULL, cluster = NULL, strata = NULL) {
  outcome_nm <- rlang::as_name(rlang::ensym(outcome))
  design_cols <- c(weights, cluster, strata)
  if (!is.null(design_cols) && (!is.character(design_cols) || any(lengths(list(weights, cluster, strata)) > 1))) {
    stop("`weights`, `cluster` and `strata` must each be a single column name (string).", call. = FALSE)
  }
  df <- drop_missing(data, c(outcome_nm, design_cols), na.rm = na.rm)
  warn_if(!success %in% unique(stats::na.omit(df[[outcome_nm]])), "The declared `success` level is not present in the data.")
  success_n <- sum(df[[outcome_nm]] == success, na.rm = TRUE)
  total_n <- sum(!is.na(df[[outcome_nm]]))
  binom <- stats::binom.test(success_n, total_n, p = p)
  prop <- suppressWarnings(stats::prop.test(success_n, total_n, p = p, correct = FALSE))
  expected_success <- total_n * p
  expected_failure <- total_n * (1 - p)
  approx_ok <- expected_success >= 5 && expected_failure >= 5
  warn_if(!approx_ok && expected_success >= 3 && expected_failure >= 3, "Proportion test approximation is borderline; exact binomial results are safer.")
  effect <- .bt(name = "Observed proportion", estimate = success_n / total_n, magnitude = NA_character_)
  plt <- if (plot) {
    counts <- c(success_n, total_n - success_n)
    pd <- .bt(level = c("success", "other"), n = counts, prop = counts / total_n)
    ggplot2::ggplot(pd, ggplot2::aes(x = .data$level, y = .data$prop, fill = .data$level)) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept = p, linetype = "dashed") +
      ggplot2::labs(title = "One-sample proportion workflow", subtitle = plot_subtitle("Exact binomial test", binom), x = NULL, y = "Proportion") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  } else NULL
  h0 <- h0_proportion(outcome_nm, p)
  primary <- if (approx_ok) prop else binom
  primary_name <- if (approx_ok) "One-sample proportion test" else "Exact binomial test"
  out <- new_testflow("proportion", "one categorical proportion", outcome_nm, data = df, descriptives = descriptives_categorical(df, outcome_nm), assumptions = assumption_checks(check_independence_note(), assumption_check("Fixed number of trials", "assumed", "The workflow treats the observed sample size as fixed."), assumption_check("Expected successes and failures", ifelse(approx_ok, "acceptable", "not acceptable"), ifelse(approx_ok, "Approximate one-sample proportion test is reasonable.", "Exact binomial test is safer for small expected counts."), details = paste0("expected_success = ", format(expected_success, digits = 3), "; expected_failure = ", format(expected_failure, digits = 3)))), recommended = list(test = primary_name), primary_test = add_null_hypothesis(safe_tidy_htest(primary, primary_name), h0), alternative_tests = list(exact_binomial = add_null_hypothesis(safe_tidy_htest(binom, "Exact binomial test"), h0), prop_test = add_null_hypothesis(safe_tidy_htest(prop, "One-sample proportion test"), h0)), effect_size = effect, plot = plt, call = match.call(), subclass = "proportion")
  if (length(design_cols)) {
    survey_res <- proportion_survey_test(df, outcome_nm, success, p, alpha, weights, cluster, strata)
    out$alternative_tests <- c(list(primary_unweighted = out$primary_test), out$alternative_tests)
    out$primary_test <- add_null_hypothesis(safe_tidy_htest(survey_res$test, "Design-based proportion test"), h0)
    out$primary_test$conf.level <- 1 - alpha
    out$recommended <- list(test = "Design-based proportion test", rationale = "Survey design columns were supplied, so inference uses the weighted estimate and design-based SE.")
    out$effect_size <- .bt(name = "Weighted proportion", estimate = survey_res$estimate, magnitude = NA_character_)
    out$survey <- survey_res$summary
    out$assumptions <- rbind(out$assumptions, assumption_check("Survey design", "assumed", "Weights, clusters and strata describe the sampling design; clusters are treated as sampled with replacement within strata.", details = paste0("design effect = ", format(survey_res$summary$deff, digits = 3), "; design df = ", survey_res$summary$df)))
  }
  out$interpretation <- make_report(out, alpha)
  out
}

proportion_survey_test <- function(df, outcome_nm, success, p, alpha, weights, cluster, strata) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop("Package `survey` is required for weighted proportions. Install it with install.packages(\"survey\").", call. = FALSE)
  }
  df$.testflow_success <- as.numeric(df[[outcome_nm]] == success)
  as_f <- function(x) if (is.null(x)) NULL else stats::as.formula(paste0("~", x))
  design <- survey::svydesign(
    ids = as_f(cluster) %||% ~1, strata = as_f(strata), weights = as_f(weights),
    data = df, nest = TRUE
  )
  est <- survey::svymean(~.testflow_success, design, deff = TRUE)
  phat <- unname(stats::coef(est))
  se <- unname(survey::SE(est))
  # Same df as survey::svyttest() (design df minus one for the mean).
  dfree <- survey::degf(design) - 1
  tstat <- (phat - p) / se
  ci <- survey::svyciprop(~.testflow_success, design, method = "logit", level = 1 - alpha)
  ci <- unname(stats::confint(ci))
  test <- structure(list(
    statistic = c(t = tstat), parameter = c(df = dfree), p.value = 2 * stats::pt(-abs(tstat), dfree),
    estimate = c(p = phat), null.value = c(p = p), conf.int = structure(as.numeric(ci), conf.level = 1 - alpha),
    alternative = "two.sided", method = "Design-based one-sample proportion test", data.name = outcome_nm
  ), class = "htest")
  deff <- unname(survey::deff(est))
  list(test = test, estimate = phat, summary = .bt(estimate = phat, se = se, conf.low = ci[1], conf.high = ci[2], df = dfree, deff = deff, n = nrow(df)))
}
