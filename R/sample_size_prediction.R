#' Sample size for developing a clinical prediction model (Riley criteria)
#'
#' Minimum sample size for a new multivariable prediction model following
#' Riley et al. (2019, 2020): the largest n meeting every criterion below.
#'
#' @param outcome `"binary"` (logistic), `"survival"` (Cox or other
#'   time-to-event), or `"continuous"` (linear regression).
#' @param parameters Number of candidate predictor parameters (degrees of
#'   freedom, counting each dummy and spline term), not number of variables.
#' @param csrsquared Anticipated Cox-Snell R-squared (binary, survival), for
#'   example from an existing model: `1 - exp(-LR / n)`.
#' @param nagrsquared Anticipated Nagelkerke R-squared (binary, survival);
#'   converted to Cox-Snell. Give one of `csrsquared` or `nagrsquared`.
#' @param rsquared Anticipated R-squared (continuous).
#' @param prevalence Outcome proportion (binary).
#' @param rate Event rate per person-time unit (survival).
#' @param timepoint Time horizon for the risk prediction (survival).
#' @param meanfup Anticipated mean follow-up per person (survival).
#' @param intercept,sd Anticipated mean outcome and its SD (continuous).
#' @param shrinkage Target expected uniform shrinkage (default 0.9, i.e. at
#'   most 10% overfitting).
#' @param mmoe Multiplicative margin of error for the mean outcome
#'   (continuous; default 1.1).
#' @return A `sample_size` object. `$criteria` holds the n, shrinkage and
#'   events (or subjects) per parameter for each criterion.
#' @details
#' Criterion 1 (shrinkage \eqn{S} of at least `shrinkage`) for binary and
#' survival outcomes, with \eqn{p} parameters and Cox-Snell \eqn{R^2_{CS}}:
#' \deqn{n_1 = \frac{p}{(S - 1)\log(1 - R^2_{CS}/S)}}
#' Criterion 2 (optimism in Nagelkerke \eqn{R^2} at most 0.05) uses the same
#' formula with \eqn{S = R^2_{CS}/(R^2_{CS} + 0.05\,\max R^2_{CS})}, where
#' \eqn{\max R^2_{CS} = 1 - \exp(2 \ln L_{null}/n)}. Criterion 3 for binary
#' outcomes estimates the overall risk within 0.05:
#' \eqn{n_3 = (1.96/0.05)^2 \phi(1-\phi)}. For survival, criterion 3 reports
#' the precision of the cumulative incidence at `timepoint`.
#' For continuous outcomes the four criteria are shrinkage, optimism in
#' \eqn{R^2}, precision of the residual SD, and precision of the mean
#' outcome (Riley et al. 2019, part I).
#' The calculations follow the published formulas and agree with
#' `pmsampsize` (tests check the binary, survival and continuous cases).
#' `pmsampsize` can also derive \eqn{R^2_{CS}} from an anticipated
#' C-statistic by simulation; that option is not provided here.
#' @references
#' Riley RD, Snell KIE, Ensor J, et al. Minimum sample size for developing a
#' multivariable prediction model: Part I - continuous outcomes. *Statistics
#' in Medicine*. 2019;38(7):1262-1275.
#'
#' Riley RD, Snell KIE, Ensor J, et al. Minimum sample size for developing a
#' multivariable prediction model: Part II - binary and time-to-event
#' outcomes. *Statistics in Medicine*. 2019;38(7):1276-1296.
#'
#' Riley RD, Ensor J, Snell KIE, et al. Calculating the sample size required
#' for developing a clinical prediction model. *BMJ*. 2020;368:m441.
#' @examples
#' samplesizeprediction("binary", parameters = 24, csrsquared = 0.288, prevalence = 0.174)
#' samplesizeprediction("survival", parameters = 25, csrsquared = 0.051,
#'                      rate = 0.065, timepoint = 2, meanfup = 2.07)
#' @export
samplesizeprediction <- function(
  outcome = c("binary", "survival", "continuous"),
  parameters,
  csrsquared = NULL,
  nagrsquared = NULL,
  rsquared = NULL,
  prevalence = NULL,
  rate = NULL,
  timepoint = NULL,
  meanfup = NULL,
  intercept = NULL,
  sd = NULL,
  shrinkage = 0.9,
  mmoe = 1.1
) {
  outcome <- match.arg(outcome)
  sample_size_validate_positive(parameters, "parameters")
  sample_size_validate_probability(shrinkage, "shrinkage")
  shrink_n <- function(r2, s) ceiling(parameters / ((s - 1) * log(1 - r2 / s)))

  if (outcome == "continuous") {
    sample_size_validate_probability(rsquared, "rsquared")
    sample_size_validate_positive(sd, "sd")
    if (is.null(intercept) || !is.finite(intercept) || intercept == 0) stop("Give a non-zero anticipated mean outcome as `intercept`.", call. = FALSE)
    return(samplesizeprediction_continuous(parameters, rsquared, intercept, sd, shrinkage, mmoe))
  }

  if (!xor(is.null(csrsquared), is.null(nagrsquared))) {
    stop("Give exactly one of `csrsquared` or `nagrsquared`.", call. = FALSE)
  }

  if (outcome == "binary") {
    sample_size_validate_probability(prevalence, "prevalence")
    max_r2 <- 1 - (prevalence^prevalence * (1 - prevalence)^(1 - prevalence))^2
    events_per_n <- prevalence
  } else {
    sample_size_validate_positive(rate, "rate")
    sample_size_validate_positive(timepoint, "timepoint")
    sample_size_validate_positive(meanfup, "meanfup")
    # Same reference cohort of 10,000 people as Riley et al. and pmsampsize.
    n_ref <- 10000
    ev <- ceiling(rate * meanfup * n_ref)
    max_r2 <- 1 - exp(2 * (ev * log(ev / n_ref) - ev) / n_ref)
    events_per_n <- rate * meanfup
  }
  r2 <- csrsquared %||% (nagrsquared * max_r2)
  sample_size_validate_probability(r2, "Cox-Snell R-squared")
  if (r2 >= max_r2) stop("The Cox-Snell R-squared exceeds its maximum (", format_sample_value(max_r2), ") for this outcome rate.", call. = FALSE)
  if (shrinkage <= r2) stop("`shrinkage` must exceed the Cox-Snell R-squared.", call. = FALSE)

  n1 <- shrink_n(r2, shrinkage)
  s2 <- r2 / (r2 + 0.05 * max_r2)
  n2 <- shrink_n(r2, s2)
  if (outcome == "binary") {
    n3 <- ceiling((1.96 / 0.05)^2 * prevalence * (1 - prevalence))
    c3_label <- "Overall risk within +/- 0.05"
    precision_note <- NULL
  } else {
    n3 <- max(n1, n2)
    pt_total <- round(meanfup * n3, 1)
    cuminc <- 1 - exp(-rate * timepoint)
    ci <- 1 - exp(-(rate + c(-1, 1) * 1.96 * sqrt(rate / pt_total)) * timepoint)
    c3_label <- "Precise overall risk at timepoint"
    precision_note <- paste0(
      "At n = ", n3, ", the cumulative incidence at t = ", format_sample_value(timepoint), " is ",
      format_sample_value(cuminc), " (95% CI ", format_sample_value(ci[1]), " to ", format_sample_value(ci[2]), ")."
    )
  }
  n_final <- max(n1, n2, n3)
  shrink_final <- max(shrinkage, s2)
  criteria <- .bt(
    criterion = c("1: shrinkage", "2: optimism in R-squared <= 0.05", paste0("3: ", c3_label), "Final"),
    n = c(n1, n2, n3, n_final),
    shrinkage = c(shrinkage, s2, shrink_final, shrink_final),
    events = c(n1, n2, n3, n_final) * events_per_n,
    epp = c(n1, n2, n3, n_final) * events_per_n / parameters
  )
  out <- sample_size_result(
    endpoint = outcome,
    design = "prediction model development",
    objective = "minimum sample size (Riley criteria)",
    method = paste0("Riley et al. (2019, 2020), ", if (outcome == "binary") "logistic model" else "time-to-event model"),
    n = n_final,
    n_adjusted = n_final,
    n_total = n_final,
    required_events = ceiling(n_final * events_per_n),
    assumptions = sample_size_assumptions(
      paste0(parameters, " candidate predictor parameters; Cox-Snell R-squared = ", format_sample_value(r2),
             " (max ", format_sample_value(max_r2), ")."),
      paste0("Criterion 1 (shrinkage >= ", format_sample_value(shrinkage), "): n = ", n1, "."),
      paste0("Criterion 2 (optimism in Nagelkerke R-squared <= 0.05, shrinkage ", format_sample_value(s2), "): n = ", n2, "."),
      paste0("Criterion 3 (", tolower(c3_label), "): n = ", n3, "."),
      precision_note
    ),
    report = paste0(
      "To develop a ", outcome, " prediction model with ", parameters, " candidate parameters and an anticipated Cox-Snell R-squared of ",
      format_sample_value(r2), ", at least ", n_final, " participants are needed (",
      format_sample_value(n_final * events_per_n), " events, ", format_sample_value(n_final * events_per_n / parameters),
      " events per parameter), following Riley et al."
    ),
    plot_data = sample_size_plot_frame(c("Criterion 1", "Criterion 2", "Criterion 3"), c(n1, n2, n3), c(n1, n2, n3))
  )
  out$criteria <- criteria
  out
}

samplesizeprediction_continuous <- function(parameters, r2, intercept, sd, shrinkage, mmoe) {
  shrink_at <- function(n) 1 + (parameters - 2) / (n * log(1 - ((r2 * (n - parameters - 1)) + parameters) / (n - 1)))
  n1 <- parameters + 2
  while (is.na(shrink_at(n1)) || shrink_at(n1) < shrinkage) n1 <- n1 + 1
  n2 <- ceiling(1 + parameters * (1 - r2) / 0.05)
  n3 <- 234 + parameters
  mean_moe <- function(n) (intercept + sqrt(sd^2 * (1 - r2) / n) * stats::qt(0.975, n - parameters - 1)) / intercept
  n4 <- max(n1, n2, n3)
  while (mean_moe(n4) > mmoe) n4 <- n4 + 1
  n_final <- max(n1, n2, n3, n4)
  criteria <- .bt(
    criterion = c("1: shrinkage", "2: optimism in R-squared <= 0.05", "3: residual SD within 10%", "4: mean outcome within margin", "Final"),
    n = c(n1, n2, n3, n4, n_final),
    shrinkage = vapply(c(n1, n2, n3, n4, n_final), shrink_at, numeric(1)),
    spp = c(n1, n2, n3, n4, n_final) / parameters
  )
  out <- sample_size_result(
    endpoint = "continuous",
    design = "prediction model development",
    objective = "minimum sample size (Riley criteria)",
    method = "Riley et al. (2019), linear regression model",
    n = n_final,
    n_adjusted = n_final,
    n_total = n_final,
    assumptions = sample_size_assumptions(
      paste0(parameters, " candidate predictor parameters; anticipated R-squared = ", format_sample_value(r2), "."),
      paste0("Criterion 1 (shrinkage >= ", format_sample_value(shrinkage), "): n = ", n1, "."),
      paste0("Criterion 2 (optimism in R-squared <= 0.05): n = ", n2, "."),
      paste0("Criterion 3 (residual SD within 10%): n = ", n3, "."),
      paste0("Criterion 4 (mean outcome within a factor ", format_sample_value(mmoe), "): n = ", n4, ".")
    ),
    report = paste0(
      "To develop a linear prediction model with ", parameters, " candidate parameters and an anticipated R-squared of ",
      format_sample_value(r2), ", at least ", n_final, " participants are needed (",
      format_sample_value(n_final / parameters), " subjects per parameter), following Riley et al."
    ),
    plot_data = sample_size_plot_frame(paste("Criterion", 1:4), c(n1, n2, n3, n4), c(n1, n2, n3, n4))
  )
  out$criteria <- criteria
  out
}
