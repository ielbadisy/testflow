#' Sample size planning for common trial designs
#'
#' @description
#' `samplesize()` dispatches to endpoint-specific planning helpers for
#' continuous, binary, survival, and ordinal planning problems.
#'
#' @param endpoint Endpoint family: `continuous`, `binary`, `survival`, or
#'   `ordinal`.
#' @param design Study design: `parallel`, `paired`, or `repeated`.
#' @param objective Planning objective: `superiority`, `noninferiority`,
#'   `equivalence`, or `precision` depending on endpoint.
#' @param ... Endpoint-specific arguments passed to the selected helper.
#' @return A `sample_size` object.
#' @references
#' Julious SA. *Sample Sizes for Clinical Trials*. Chapman & Hall/CRC; 2010.
#' @export
samplesize <- function(endpoint = c("continuous", "binary", "survival", "ordinal"), design = c("parallel", "paired", "repeated"), objective = c("superiority", "noninferiority", "equivalence"), ...) {
  endpoint <- match.arg(endpoint)
  design <- match.arg(design)
  objective <- match.arg(objective)
  switch(
    endpoint,
    continuous = samplesizecontinuous(design = design, objective = objective, ...),
    binary = samplesizebinary(design = design, objective = objective, ...),
    survival = samplesizesurvival(design = design, objective = objective, ...),
    ordinal = samplesizeordinal(design = design, objective = objective, ...)
  )
}

#' Sample size for continuous endpoints
#'
#' @param design `parallel`, `paired`, or `repeated`.
#' @param objective `superiority`, `noninferiority`, or `equivalence`.
#' @param delta Effect size for superiority or equivalence margin.
#' @param sd Common between-subject standard deviation.
#' @param sd_diff Standard deviation of paired differences. If omitted for
#'   paired designs, `sd` is used.
#' @param expected_difference Planned treatment difference.
#' @param margin Non-inferiority or equivalence margin.
#' @param alpha Type I error rate.
#' @param power Target power.
#' @param allocation Allocation ratio `n_B / n_A` for parallel designs.
#' @param dropout Expected dropout proportion.
#' @param n_time Number of repeated measures.
#' @param correlation Within-subject correlation for repeated measures when
#'   `design = "repeated"` and `n_time > 2`.
#' @param method `"normal"` (default) uses the normal approximation;
#'   `"t"` returns the smallest n whose noncentral-t power reaches `power`
#'   (as `stats::power.t.test()`).
#' @param baseline_cor Correlation between the baseline and follow-up
#'   measurement for an ANCOVA analysis of a parallel design. The planning
#'   SD becomes \eqn{\sigma\sqrt{1 - \rho^2}} (Borm et al. 2007). Default
#'   `0` (no baseline adjustment, results unchanged).
#' @return A `sample_size` object.
#' @details
#' With \eqn{r} the allocation ratio \eqn{n_B / n_A}, \eqn{z_{1-\beta}} and
#' \eqn{z_{1-\alpha/2}} (or \eqn{z_{1-\alpha}} for one-sided designs) standard
#' normal quantiles:
#'
#' Parallel superiority (equal allocation):
#' \deqn{n_{per\ group} = \frac{2\sigma^2(z_{1-\beta}+z_{1-\alpha/2})^2}{\Delta^2}}
#'
#' Parallel non-inferiority, with \eqn{d_{NI} = \Delta + \delta}:
#' \deqn{n_{per\ group} = \frac{2\sigma^2(z_{1-\beta}+z_{1-\alpha})^2}{d_{NI}^2}}
#'
#' Parallel equivalence, with \eqn{d_{EQ} = \delta - |\Delta|}:
#' \deqn{n_{per\ group} \approx \frac{2\sigma^2(z_{1-\beta}+z_{1-\alpha})^2}{d_{EQ}^2}}
#'
#' Paired designs use the same forms with the standard deviation of the
#' paired differences (`sd_diff`) in place of \eqn{\sigma}. Repeated designs
#' (`n_time > 2`) plan on an effective standard deviation derived from a
#' compound-symmetry working correlation.
#' @references
#' Julious SA. *Sample Sizes for Clinical Trials*. Chapman & Hall/CRC; 2010.
#' @export
samplesizecontinuous <- function(
  design = c("parallel", "paired", "repeated"),
  objective = c("superiority", "noninferiority", "equivalence"),
  delta = NULL,
  sd = NULL,
  sd_diff = NULL,
  expected_difference = 0,
  margin = NULL,
  alpha = 0.05,
  power = 0.90,
  allocation = 1,
  dropout = 0,
  n_time = 2,
  correlation = 0.5,
  method = c("normal", "t"),
  baseline_cor = 0
) {
  design <- match.arg(design)
  objective <- match.arg(objective)
  method <- match.arg(method)
  if (baseline_cor != 0) {
    if (design != "parallel") stop("`baseline_cor` applies to parallel designs only.", call. = FALSE)
    if (abs(baseline_cor) >= 1) stop("`baseline_cor` must be in (-1, 1).", call. = FALSE)
    sample_size_validate_positive(sd, "sd")
    sd_a <- sd * sqrt(1 - baseline_cor^2)
    res <- samplesizecontinuous(design, objective, delta, sd_a, sd_diff, expected_difference, margin, alpha, power,
                                allocation, dropout, n_time, correlation, method)
    res$method <- paste0(res$method, ", ANCOVA on baseline"); res$inputs <- ssinputs(samplesizecontinuous, environment())
    res$assumptions <- rbind(res$assumptions, sample_size_assumptions(paste0(
      "ANCOVA: baseline-outcome correlation = ", format_sample_value(baseline_cor), ", SD ", format_sample_value(sd),
      " reduced to ", format_sample_value(sd_a), "."))[1, ])
    res$assumptions$name[nrow(res$assumptions)] <- paste0("Planning note ", nrow(res$assumptions))
    res$report <- paste0(res$report, " The SD was reduced from ", format_sample_value(sd), " to ", format_sample_value(sd_a),
                         " for ANCOVA with baseline correlation ", format_sample_value(baseline_cor), ".")
    return(res)
  }
  use_t <- identical(method, "t")
  t_label <- function(x) if (use_t) paste0(sub(" normal approximation$| approximation$", "", x), " (noncentral t)") else x
  sample_size_validate_probability(alpha, "alpha")
  sample_size_validate_probability(power, "power")
  sample_size_validate_dropout(dropout)
  sample_size_validate_positive(allocation, "allocation")
  sample_size_validate_probability(correlation, "correlation", allow_zero = TRUE)
  if (design == "repeated" && n_time < 2) {
    stop("Repeated-measures sample size requires at least 2 time points.", call. = FALSE)
  }

  z_power <- stats::qnorm(power)
  z_alpha_two <- stats::qnorm(1 - alpha / 2)
  z_alpha_one <- stats::qnorm(1 - alpha)
  paired_design <- design %in% c("paired", "repeated")
  repeated_design <- identical(design, "repeated") && n_time > 2
  design_label <- if (repeated_design) {
    paste0("repeated (", n_time, " time points)")
  } else if (design == "repeated") {
    "paired (two time points)"
  } else {
    "paired"
  }
  planning_label <- if (repeated_design) "subjects" else "pairs"

  if (paired_design) {
    scale_sd <- sd_diff %||% sd
    sample_size_validate_positive(scale_sd, "sd_diff")
    if (repeated_design) {
      scale_sd <- sample_size_repeated_effective_sd(scale_sd, n_time = n_time, correlation = correlation)
    }

    if (objective == "superiority") {
      sample_size_validate_positive(delta, "delta")
      n_raw <- ((z_alpha_two + z_power)^2 * scale_sd^2) / delta^2
      if (use_t) n_raw <- sample_size_continuous_t(n_raw, "superiority", delta, scale_sd, alpha, power, paired = TRUE)
      n_main <- sample_size_round(n_raw)
      n_adj <- sample_size_apply_dropout(n_main, dropout)
      return(sample_size_result(
        endpoint = "continuous",
        design = design_label,
        objective = objective,
        method = t_label(if (repeated_design) "repeated-measures normal approximation" else "paired normal approximation"),
        n = n_raw,
        n_adjusted = n_adj,
        n_total = n_adj,
        assumptions = sample_size_assumptions(
          if (repeated_design) {
            paste0("Repeated measurements from the same subject are assumed (", n_time, " time points).")
          } else {
            "Paired observations are assumed."
          },
          paste0("Effective SD = ", format_sample_value(scale_sd), "."),
          if (repeated_design) paste0("Within-subject correlation = ", format_sample_value(correlation), ".") else NULL
        ),
        report = paste0(
          if (repeated_design) "Repeated-measures continuous superiority sample size was calculated with delta = " else "Paired continuous superiority sample size was calculated with delta = ",
          format_sample_value(delta), ", sd_diff = ", format_sample_value(scale_sd),
          ", alpha = ", format_sample_value(alpha), ", power = ", format_sample_value(power),
          "; required ", planning_label, " = ", format_sample_count(n_adj), "."
        ),
        plot_data = sample_size_plot_frame(paste0("Required ", planning_label), n_raw, n_adj),
        curve_data = sample_size_curve_bundle(
          n_target = n_adj,
          power_target = power,
          n_grid = sample_size_curve_grid(n_adj),
          power_grid = sample_size_curve_power_continuous_paired_superiority(sample_size_curve_grid(n_adj), delta, scale_sd, alpha = alpha),
          unit_label = "sample size",
          power_label = "test power = 1 - beta"
        )
      ))
    }

    if (objective == "noninferiority") {
      sample_size_validate_positive(delta, "margin")
      d_ni <- expected_difference + delta
      if (d_ni <= 0) {
        stop("The non-inferiority distance must be positive: expected_difference + margin > 0.", call. = FALSE)
      }
      n_raw <- ((z_power + z_alpha_one)^2 * scale_sd^2) / d_ni^2
      if (use_t) n_raw <- sample_size_continuous_t(n_raw, "noninferiority", d_ni, scale_sd, alpha, power, paired = TRUE)
      n_main <- sample_size_round(n_raw)
      n_adj <- sample_size_apply_dropout(n_main, dropout)
      return(sample_size_result(
        endpoint = "continuous",
        design = design_label,
        objective = objective,
        method = t_label(if (repeated_design) "repeated-measures non-inferiority normal approximation" else "paired non-inferiority normal approximation"),
        n = n_raw,
        n_adjusted = n_adj,
        n_total = n_adj,
        assumptions = sample_size_assumptions(
          if (repeated_design) {
            paste0("Repeated measurements from the same subject are assumed (", n_time, " time points).")
          } else {
            "Paired observations are assumed."
          },
          paste0("Effective SD = ", format_sample_value(scale_sd), "."),
          if (repeated_design) paste0("Within-subject correlation = ", format_sample_value(correlation), ".") else NULL,
          paste0("Distance to the non-inferiority boundary = ", format_sample_value(d_ni), ".")
        ),
        report = paste0(
          if (repeated_design) "Repeated-measures continuous non-inferiority sample size was calculated with expected_difference = " else "Paired continuous non-inferiority sample size was calculated with expected_difference = ",
          format_sample_value(expected_difference), ", margin = ", format_sample_value(delta),
          ", sd_diff = ", format_sample_value(scale_sd), "; required ", planning_label, " = ", format_sample_count(n_adj), "."
        ),
        plot_data = sample_size_plot_frame(paste0("Required ", planning_label), n_raw, n_adj)
      ))
    }

    if (objective == "equivalence") {
      sample_size_validate_positive(delta, "margin")
      d_eq <- delta - abs(expected_difference)
      if (d_eq <= 0) {
        stop("The equivalence distance must be positive: margin must exceed the absolute expected difference.", call. = FALSE)
      }
      if (isTRUE(all.equal(expected_difference, 0))) {
        z_power_eq <- stats::qnorm(1 - (1 - power) / 2)
        n_raw <- 2 * scale_sd^2 * (z_power_eq + z_alpha_one)^2 / delta^2
      } else {
        n_raw <- ((z_power + z_alpha_one)^2 * scale_sd^2) / d_eq^2
      }
      if (use_t) n_raw <- sample_size_continuous_t(n_raw, "equivalence", expected_difference, scale_sd, alpha, power, paired = TRUE, margin = delta)
      n_main <- sample_size_round(n_raw)
      n_adj <- sample_size_apply_dropout(n_main, dropout)
      return(sample_size_result(
        endpoint = "continuous",
        design = design_label,
        objective = objective,
        method = t_label(if (repeated_design) "repeated-measures equivalence approximation" else "paired equivalence approximation"),
        n = n_raw,
        n_adjusted = n_adj,
        n_total = n_adj,
        assumptions = sample_size_assumptions(
          if (repeated_design) {
            paste0("Repeated measurements from the same subject are assumed (", n_time, " time points).")
          } else {
            "Paired observations are assumed."
          },
          paste0("Effective SD = ", format_sample_value(scale_sd), "."),
          if (repeated_design) paste0("Within-subject correlation = ", format_sample_value(correlation), ".") else NULL,
          paste0("Distance to the nearest equivalence boundary = ", format_sample_value(d_eq), ".")
        ),
        report = paste0(
          if (repeated_design) "Repeated-measures continuous equivalence sample size was calculated with expected_difference = " else "Paired continuous equivalence sample size was calculated with expected_difference = ",
          format_sample_value(expected_difference), ", margin = ", format_sample_value(delta),
          ", sd_diff = ", format_sample_value(scale_sd), "; required ", planning_label, " = ", format_sample_count(n_adj), "."
        ),
        plot_data = sample_size_plot_frame(paste0("Required ", planning_label), n_raw, n_adj)
      ))
    }

    stop("Unsupported objective for paired continuous sample size.", call. = FALSE)
  }

  sample_size_validate_positive(sd, "sd")

  if (objective == "superiority") {
    sample_size_validate_positive(delta, "delta")
    n_a_raw <- ((allocation + 1) * sd^2 * (z_alpha_two + z_power)^2) / (allocation * delta^2)
    if (use_t) n_a_raw <- sample_size_continuous_t(n_a_raw, "superiority", delta, sd, alpha, power, allocation = allocation)
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "continuous",
      design = "parallel",
      objective = objective,
      method = t_label("parallel normal approximation"),
      n = n_a_raw,
      n_adjusted = c(A = n_a_adj, B = n_b_adj),
      n_per_group = c(A = n_a_adj, B = n_b_adj),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        paste0("Allocation ratio B/A = ", format_sample_value(allocation), ".")
      ),
      report = paste0(
        "Parallel continuous superiority sample size was calculated with delta = ",
        format_sample_value(delta), ", sd = ", format_sample_value(sd), ", allocation ratio = ",
        format_sample_value(allocation), "; required per-group sizes = A ", format_sample_count(n_a_adj),
        ", B ", format_sample_count(n_b_adj), "."
      ),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_a_raw, n_b_raw), c(n_a_adj, n_b_adj))
    ))
  }

  if (objective == "noninferiority") {
    sample_size_validate_positive(delta, "margin")
    d_ni <- expected_difference + delta
    if (d_ni <= 0) {
      stop("The non-inferiority distance must be positive: expected_difference + margin > 0.", call. = FALSE)
    }
    n_a_raw <- ((allocation + 1) * sd^2 * (z_power + z_alpha_one)^2) / (allocation * d_ni^2)
    if (use_t) n_a_raw <- sample_size_continuous_t(n_a_raw, "noninferiority", d_ni, sd, alpha, power, allocation = allocation)
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "continuous",
      design = "parallel",
      objective = objective,
      method = t_label("parallel non-inferiority normal approximation"),
      n = n_a_raw,
      n_adjusted = c(A = n_a_adj, B = n_b_adj),
      n_per_group = c(A = n_a_adj, B = n_b_adj),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        paste0("Distance to the non-inferiority boundary = ", format_sample_value(d_ni), ".")
      ),
      report = paste0(
        "Parallel continuous non-inferiority sample size was calculated with expected_difference = ",
        format_sample_value(expected_difference), ", margin = ", format_sample_value(delta),
        ", sd = ", format_sample_value(sd), "; required per-group sizes = A ", format_sample_count(n_a_adj),
        ", B ", format_sample_count(n_b_adj), "."
      ),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_a_raw, n_b_raw), c(n_a_adj, n_b_adj))
    ))
  }

  if (objective == "equivalence") {
    sample_size_validate_positive(delta, "margin")
    d_eq <- delta - abs(expected_difference)
    if (d_eq <= 0) {
      stop("The equivalence distance must be positive: margin must exceed the absolute expected difference.", call. = FALSE)
    }
    if (isTRUE(all.equal(expected_difference, 0))) {
      z_power_eq <- stats::qnorm(1 - (1 - power) / 2)
      n_a_raw <- ((allocation + 1) * sd^2 * (z_power_eq + z_alpha_one)^2) / (allocation * delta^2)
    } else {
      n_a_raw <- ((allocation + 1) * sd^2 * (z_power + z_alpha_one)^2) / (allocation * d_eq^2)
    }
    if (use_t) n_a_raw <- sample_size_continuous_t(n_a_raw, "equivalence", expected_difference, sd, alpha, power, allocation = allocation, margin = delta)
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "continuous",
      design = "parallel",
      objective = objective,
      method = t_label("parallel equivalence approximation"),
      n = n_a_raw,
      n_adjusted = c(A = n_a_adj, B = n_b_adj),
      n_per_group = c(A = n_a_adj, B = n_b_adj),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        paste0("Distance to the nearest equivalence boundary = ", format_sample_value(d_eq), ".")
      ),
      report = paste0(
        "Parallel continuous equivalence sample size was calculated with expected_difference = ",
        format_sample_value(expected_difference), ", margin = ", format_sample_value(delta),
        ", sd = ", format_sample_value(sd), "; required per-group sizes = A ", format_sample_count(n_a_adj),
        ", B ", format_sample_count(n_b_adj), "."
      ),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_a_raw, n_b_raw), c(n_a_adj, n_b_adj))
    ))
  }

  stop("Unsupported objective for continuous sample size.", call. = FALSE)
}

#' Sample size for binary endpoints
#'
#' @param design `parallel`, `paired`, or `repeated`.
#' @param objective `superiority`, `noninferiority`, or `equivalence`.
#' @param p1 Probability in arm A or first paired state.
#' @param p2 Probability in arm B or second paired state.
#' @param margin Non-inferiority or equivalence margin on the risk-difference scale.
#' @param method Variance estimator for parallel superiority planning:
#'   `pooled` uses the null (pooled-proportion) variance and `anticipated`
#'   uses the alternative-hypothesis variance. Ignored for other objectives.
#' @param discordant_or Discordant odds ratio for paired binary superiority.
#' @param discordance_rate Proportion of discordant pairs.
#' @param p10 Probability of A-only response for paired binary designs.
#' @param p01 Probability of B-only response for paired binary designs.
#' @param alpha Type I error rate.
#' @param power Target power.
#' @param allocation Allocation ratio `n_B / n_A` for parallel designs.
#' @param dropout Expected dropout proportion.
#' @param n_time Number of repeated measures.
#' @return A `sample_size` object.
#' @details
#' With \eqn{\bar p = (p_A + r p_B)/(1+r)} the pooled planning proportion and
#' \eqn{r} the allocation ratio:
#'
#' Parallel superiority, pooled variance:
#' \deqn{n_A = \frac{\left[z_{1-\alpha/2}\sqrt{(1+1/r)\bar p(1-\bar p)}+z_{1-\beta}\sqrt{p_A(1-p_A)+p_B(1-p_B)/r}\right]^2}{(p_A-p_B)^2}}
#'
#' Parallel superiority, anticipated-response variance:
#' \deqn{n_A \approx \frac{(z_{1-\alpha/2}+z_{1-\beta})^2\left[p_A(1-p_A)+p_B(1-p_B)/r\right]}{(p_A-p_B)^2}}
#'
#' Non-inferiority, with \eqn{d_{NI} = (p_A-p_B)+\delta}:
#' \deqn{n_A \approx \frac{(z_{1-\alpha}+z_{1-\beta})^2\left[p_A(1-p_A)+p_B(1-p_B)/r\right]}{d_{NI}^2}}
#'
#' Equivalence, with \eqn{d_{EQ} = \delta - |p_A-p_B|}:
#' \deqn{n_A \approx \frac{(z_{1-\alpha}+z_{1-\beta})^2\left[p_A(1-p_A)+p_B(1-p_B)/r\right]}{d_{EQ}^2}}
#'
#' Equivalence special case \eqn{p_A = p_B = p}:
#' \deqn{n_A = \frac{(r+1)(z_{1-\beta/2}+z_{1-\alpha})^2p(1-p)}{r\delta^2}}
#'
#' Paired (discordant-pairs) superiority, with discordant odds ratio
#' \eqn{OR_D = p_{10}/p_{01}} and discordance rate \eqn{\lambda_D = p_{10}+p_{01}}:
#' \deqn{n_{total} = \left\lceil\frac{(z_{1-\alpha/2}+z_{1-\beta})^2(OR_D+1)^2}{(OR_D-1)^2\lambda_D}\right\rceil}
#' @references
#' Julious SA. *Sample Sizes for Clinical Trials*. Chapman & Hall/CRC; 2010.
#' @export
samplesizebinary <- function(
  design = c("parallel", "paired", "repeated"),
  objective = c("superiority", "noninferiority", "equivalence"),
  p1,
  p2,
  margin = NULL,
  method = c("pooled", "anticipated"),
  discordant_or = NULL,
  discordance_rate = NULL,
  p10 = NULL,
  p01 = NULL,
  alpha = 0.05,
  power = 0.90,
  allocation = 1,
  dropout = 0,
  n_time = 2
) {
  design <- match.arg(design)
  objective <- match.arg(objective)
  method <- match.arg(method)
  sample_size_validate_probability(alpha, "alpha")
  sample_size_validate_probability(power, "power")
  sample_size_validate_dropout(dropout)
  sample_size_validate_positive(allocation, "allocation")
  if (design == "repeated" && n_time > 2) {
    stop("Repeated-measures binary sample size is currently supported only when `n_time = 2`.", call. = FALSE)
  }

  z_power <- stats::qnorm(power)
  z_alpha_two <- stats::qnorm(1 - alpha / 2)
  z_alpha_one <- stats::qnorm(1 - alpha)

  if (design %in% c("paired", "repeated")) {
    if (!is.null(p10) || !is.null(p01)) {
      if (is.null(p10) || is.null(p01)) {
        stop("Provide both `p10` and `p01` for paired binary planning.", call. = FALSE)
      }
      sample_size_validate_probability(p10, "p10")
      sample_size_validate_probability(p01, "p01")
      discordant_or <- p10 / p01
      discordance_rate <- p10 + p01
    }
    sample_size_validate_positive(discordant_or, "discordant_or")
    sample_size_validate_probability(discordance_rate, "discordance_rate")
    if (isTRUE(all.equal(discordant_or, 1))) {
      stop("The discordant odds ratio must differ from 1.", call. = FALSE)
    }
    n_discordant_raw <- ((z_alpha_two + z_power)^2 * (discordant_or + 1)^2) / (discordant_or - 1)^2
    n_total_raw <- n_discordant_raw / discordance_rate
    n_main <- sample_size_round(n_total_raw)
    n_adj <- sample_size_apply_dropout(n_main, dropout)
    return(sample_size_result(
      endpoint = "binary",
      design = if (design == "repeated") "paired (two time points)" else "paired",
      objective = "superiority",
      method = "discordant pairs approximation",
      n = n_total_raw,
      n_adjusted = n_adj,
      n_total = n_adj,
      assumptions = sample_size_assumptions(
        "Matched pairs are assumed.",
        paste0("Discordant odds ratio = ", format_sample_value(discordant_or), "."),
        paste0("Discordance rate = ", format_sample_value(discordance_rate), ".")
      ),
      report = paste0(
        "Paired binary superiority sample size was calculated with discordant_or = ",
        format_sample_value(discordant_or), ", discordance_rate = ", format_sample_value(discordance_rate),
        "; required pairs = ", format_sample_count(n_adj), "."
      ),
      plot_data = sample_size_plot_frame("Required pairs", n_total_raw, n_adj)
    ))
  }

  sample_size_validate_probability(p1, "p1")
  sample_size_validate_probability(p2, "p2")
  delta <- p1 - p2

  if (objective == "superiority") {
    if (isTRUE(all.equal(delta, 0))) {
      stop("The superiority risk difference must differ from 0.", call. = FALSE)
    }
    pooled_p <- (p1 + allocation * p2) / (1 + allocation)
    alt_var <- p1 * (1 - p1) + p2 * (1 - p2) / allocation
    n_pooled <- ((z_alpha_two * sqrt((1 + 1 / allocation) * pooled_p * (1 - pooled_p)) + z_power * sqrt(alt_var))^2) / delta^2
    n_alt <- ((z_alpha_two + z_power)^2 * alt_var) / delta^2
    n_a_raw <- if (method == "pooled") n_pooled else n_alt
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "binary",
      design = "parallel",
      objective = objective,
      method = if (method == "pooled") "parallel pooled normal approximation" else "parallel anticipated-response normal approximation",
      n = n_a_raw,
      n_adjusted = c(A = n_a_adj, B = n_b_adj),
      n_per_group = c(A = n_a_adj, B = n_b_adj),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        paste0("Risk difference = ", format_sample_value(delta), ".")
      ),
      report = paste0(
        "Parallel binary superiority sample size was calculated with p1 = ", format_sample_value(p1),
        ", p2 = ", format_sample_value(p2), ", allocation ratio = ", format_sample_value(allocation),
        "; required per-group sizes = A ", format_sample_count(n_a_adj), ", B ", format_sample_count(n_b_adj), "."
      ),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_a_raw, n_b_raw), c(n_a_adj, n_b_adj))
    ))
  }

  if (objective == "noninferiority") {
    sample_size_validate_positive(margin, "margin")
    d_ni <- delta + margin
    if (d_ni <= 0) {
      stop("The non-inferiority distance must be positive: (p1 - p2) + margin > 0.", call. = FALSE)
    }
    n_a_raw <- ((z_power + z_alpha_one)^2 * (p1 * (1 - p1) + p2 * (1 - p2) / allocation)) / d_ni^2
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "binary",
      design = "parallel",
      objective = objective,
      method = "parallel non-inferiority normal approximation",
      n = n_a_raw,
      n_adjusted = c(A = n_a_adj, B = n_b_adj),
      n_per_group = c(A = n_a_adj, B = n_b_adj),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        paste0("Distance to the non-inferiority boundary = ", format_sample_value(d_ni), ".")
      ),
      report = paste0(
        "Parallel binary non-inferiority sample size was calculated with p1 = ", format_sample_value(p1),
        ", p2 = ", format_sample_value(p2), ", margin = ", format_sample_value(margin),
        "; required per-group sizes = A ", format_sample_count(n_a_adj), ", B ", format_sample_count(n_b_adj), "."
      ),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_a_raw, n_b_raw), c(n_a_adj, n_b_adj))
    ))
  }

  if (objective == "equivalence") {
    sample_size_validate_positive(margin, "margin")
    d_eq <- margin - abs(delta)
    if (d_eq <= 0) {
      stop("The equivalence distance must be positive: margin must exceed the absolute risk difference.", call. = FALSE)
    }
    if (isTRUE(all.equal(delta, 0))) {
      z_power_eq <- stats::qnorm(1 - (1 - power) / 2)
      n_a_raw <- ((allocation + 1) * (z_power_eq + z_alpha_one)^2 * p1 * (1 - p1)) / (allocation * margin^2)
    } else {
      n_a_raw <- ((z_power + z_alpha_one)^2 * (p1 * (1 - p1) + p2 * (1 - p2) / allocation)) / d_eq^2
    }
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "binary",
      design = "parallel",
      objective = objective,
      method = "parallel equivalence normal approximation",
      n = n_a_raw,
      n_adjusted = c(A = n_a_adj, B = n_b_adj),
      n_per_group = c(A = n_a_adj, B = n_b_adj),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        paste0("Distance to the nearest equivalence boundary = ", format_sample_value(d_eq), ".")
      ),
      report = paste0(
        "Parallel binary equivalence sample size was calculated with p1 = ", format_sample_value(p1),
        ", p2 = ", format_sample_value(p2), ", margin = ", format_sample_value(margin),
        "; required per-group sizes = A ", format_sample_count(n_a_adj), ", B ", format_sample_count(n_b_adj), "."
      ),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_a_raw, n_b_raw), c(n_a_adj, n_b_adj))
    ))
  }

  stop("Unsupported objective for binary sample size.", call. = FALSE)
}

#' Sample size for survival endpoints
#'
#' @param design `parallel` only.
#' @param objective `superiority`, `noninferiority`, or `equivalence`.
#' @param hr Hazard ratio.
#' @param margin_hr Non-inferiority hazard-ratio margin.
#' @param lower Lower equivalence bound on the hazard-ratio scale.
#' @param upper Upper equivalence bound on the hazard-ratio scale.
#' @param survival_a Survival probability in arm A at the planning horizon.
#' @param survival_b Survival probability in arm B at the planning horizon.
#'   Default `NULL` uses `survival_a^hr` (proportional hazards) when `survival_a`
#'   is given; with neither, the size is the number of events.
#' @param alpha Type I error rate.
#' @param power Target power.
#' @param allocation Allocation ratio `n_B / n_A`.
#' @param dropout Expected dropout proportion.
#' @param method `exponential` or `ph_only`.
#' @param accrual_duration Uniform accrual (enrollment) duration. When
#'   supplied together with `follow_up`, event probabilities account for
#'   staggered enrollment instead of assuming every subject is followed for
#'   the full study duration. Requires `survival_a`/`survival_b`.
#' @param follow_up Additional follow-up duration after accrual ends. The
#'   total study duration is `accrual_duration + follow_up`.
#' @return A `sample_size` object.
#' @details
#' With \eqn{z_{1-\beta}} and \eqn{z_{1-\alpha/2}} (or \eqn{z_{1-\alpha}} for
#' one-sided designs) standard normal quantiles and \eqn{HR} the planned
#' hazard ratio:
#'
#' Superiority, exponential survival:
#' \deqn{D_{total} = \frac{4(z_{1-\alpha/2}+z_{1-\beta})^2}{[\log(HR)]^2}}
#'
#' Superiority, proportional-hazards-only approximation:
#' \deqn{D_{total} = \frac{2(HR+1)^2(z_{1-\alpha/2}+z_{1-\beta})^2}{(HR-1)^2}}
#'
#' Non-inferiority, with \eqn{d_{NI} = \log(HR_M) - \log(HR)}:
#' \deqn{D_{total} \approx \frac{4(z_{1-\alpha}+z_{1-\beta})^2}{d_{NI}^2}}
#'
#' Equivalence, with log-hazard-ratio bounds \eqn{\theta_L, \theta_U} and
#' \eqn{d_{EQ} = \min(\theta-\theta_L,\ \theta_U-\theta)}:
#' \deqn{D_{total} \approx \frac{4(z_{1-\alpha}+z_{1-\beta})^2}{d_{EQ}^2}}
#'
#' Required total events \eqn{D_{total}} are converted to a total sample size
#' using the planned survival probabilities \eqn{S_A(t), S_B(t)} under equal
#' allocation:
#' \deqn{N_{total} = \frac{2D_{total}}{(1-S_A(t))+(1-S_B(t))}}
#'
#' When `accrual_duration` (\eqn{R}) and `follow_up` are supplied, subjects
#' are assumed to enroll uniformly over \eqn{[0, R]} and survival is assumed
#' exponential with hazard implied by \eqn{S_A(t)}/\eqn{S_B(t)} at the total
#' study duration \eqn{T = R + \text{follow\_up}}. The event probability for
#' each arm becomes the accrual-averaged
#' \deqn{\bar q_g = 1 - \frac{e^{-\lambda_g(T-R)} - e^{-\lambda_g T}}{\lambda_g R}}
#' in place of the flat \eqn{1-S_g(t)} above; this reduces to \eqn{1-S_g(T)}
#' as \eqn{R \to 0} (instantaneous accrual).
#' @references
#' Julious SA. *Sample Sizes for Clinical Trials*. Chapman & Hall/CRC; 2010.
#' @export
samplesizesurvival <- function(
  design = c("parallel"),
  objective = c("superiority", "noninferiority", "equivalence"),
  hr,
  margin_hr = NULL,
  lower = NULL,
  upper = NULL,
  survival_a = NULL,
  survival_b = NULL,
  alpha = 0.05,
  power = 0.90,
  allocation = 1,
  dropout = 0,
  method = c("exponential", "ph_only"),
  accrual_duration = NULL,
  follow_up = NULL
) {
  design <- match.arg(design)
  objective <- match.arg(objective)
  method <- match.arg(method)
  if (!identical(design, "parallel")) {
    stop("Only parallel survival sample size is implemented.", call. = FALSE)
  }
  sample_size_validate_probability(alpha, "alpha")
  sample_size_validate_probability(power, "power")
  sample_size_validate_dropout(dropout)
  sample_size_validate_positive(allocation, "allocation")
  if (!isTRUE(all.equal(allocation, 1))) {
    stop("Survival sample size currently assumes equal allocation.", call. = FALSE)
  }
  sample_size_validate_positive(hr, "hr")
  # Without survival_b the experimental arm follows proportional hazards; before, n silently equalled the events.
  if (!is.null(survival_a) && is.null(survival_b)) survival_b <- survival_a^hr

  z_power <- stats::qnorm(power)
  z_alpha_two <- stats::qnorm(1 - alpha / 2)
  z_alpha_one <- stats::qnorm(1 - alpha)
  accrual_note <- if (!is.null(accrual_duration) && !is.null(follow_up)) {
    paste0("Uniform accrual over ", format_sample_value(accrual_duration), " time units, plus ", format_sample_value(follow_up), " additional follow-up, is assumed.")
  } else {
    NULL
  }

  if (objective == "superiority") {
    if (isTRUE(all.equal(log(hr), 0))) {
      stop("The hazard ratio must differ from 1 for superiority planning.", call. = FALSE)
    }
    if (method == "exponential") {
      events_per_arm <- 2 * (z_alpha_two + z_power)^2 / (log(hr)^2)
    } else {
      events_per_arm <- ((hr + 1)^2 * (z_alpha_two + z_power)^2) / (2 * (hr - 1)^2)
    }
    total_events <- 2 * events_per_arm
    n_total_raw <- sample_size_events_to_n(total_events, survival_a, survival_b, accrual_duration, follow_up)
    n_total <- sample_size_round(n_total_raw)
    n_adj <- sample_size_apply_dropout(n_total, dropout)
    return(sample_size_result(
      endpoint = "survival",
      design = "parallel",
      objective = objective,
      method = paste0("survival ", method, " event approximation"),
      n = n_total_raw,
      n_adjusted = n_adj,
      n_total = n_adj,
      required_events = sample_size_round(total_events),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        "Event probabilities are derived from the planned survival probabilities.",
        accrual_note
      ),
      report = paste0(
        "Survival superiority sample size was calculated with hr = ", format_sample_value(hr),
        "; required total events = ", format_sample_count(sample_size_round(total_events)),
        ", estimated total sample size = ", format_sample_count(n_adj), "."
      ),
      plot_data = sample_size_plot_frame("Required total sample size", n_total_raw, n_adj)
    ))
  }

  if (objective == "noninferiority") {
    sample_size_validate_positive(margin_hr, "margin_hr")
    d_ni <- log(margin_hr) - log(hr)
    if (d_ni <= 0) {
      stop("The non-inferiority distance must be positive: log(margin_hr) - log(hr) > 0.", call. = FALSE)
    }
    total_events <- 4 * (z_alpha_one + z_power)^2 / (d_ni^2)
    n_total_raw <- sample_size_events_to_n(total_events, survival_a, survival_b, accrual_duration, follow_up)
    n_total <- sample_size_round(n_total_raw)
    n_adj <- sample_size_apply_dropout(n_total, dropout)
    return(sample_size_result(
      endpoint = "survival",
      design = "parallel",
      objective = objective,
      method = "survival non-inferiority event approximation",
      n = n_total_raw,
      n_adjusted = n_adj,
      n_total = n_adj,
      required_events = sample_size_round(total_events),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        paste0("Distance to the non-inferiority boundary = ", format_sample_value(d_ni), "."),
        accrual_note
      ),
      report = paste0(
        "Survival non-inferiority sample size was calculated with hr = ", format_sample_value(hr),
        ", margin_hr = ", format_sample_value(margin_hr), "; required total events = ",
        format_sample_count(sample_size_round(total_events)), ", estimated total sample size = ",
        format_sample_count(n_adj), "."
      ),
      plot_data = sample_size_plot_frame("Required total sample size", n_total_raw, n_adj)
    ))
  }

  if (objective == "equivalence") {
    if (is.null(lower) || is.null(upper)) {
      stop("Both `lower` and `upper` bounds are required.", call. = FALSE)
    }
    sample_size_validate_positive(lower, "lower")
    sample_size_validate_positive(upper, "upper")
    theta <- log(hr)
    d_eq <- min(theta - log(lower), log(upper) - theta)
    if (d_eq <= 0) {
      stop("The equivalence distance must be positive: the planned hazard ratio must lie inside the bounds.", call. = FALSE)
    }
    total_events <- 4 * (z_alpha_one + z_power)^2 / (d_eq^2)
    n_total_raw <- sample_size_events_to_n(total_events, survival_a, survival_b, accrual_duration, follow_up)
    n_total <- sample_size_round(n_total_raw)
    n_adj <- sample_size_apply_dropout(n_total, dropout)
    return(sample_size_result(
      endpoint = "survival",
      design = "parallel",
      objective = objective,
      method = "survival equivalence event approximation",
      n = n_total_raw,
      n_adjusted = n_adj,
      n_total = n_adj,
      required_events = sample_size_round(total_events),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        paste0("Distance to the nearest equivalence boundary = ", format_sample_value(d_eq), "."),
        accrual_note
      ),
      report = paste0(
        "Survival equivalence sample size was calculated with hr = ", format_sample_value(hr),
        "; required total events = ", format_sample_count(sample_size_round(total_events)),
        ", estimated total sample size = ", format_sample_count(n_adj), "."
      ),
      plot_data = sample_size_plot_frame("Required total sample size", n_total_raw, n_adj)
    ))
  }

  stop("Unsupported objective for survival sample size.", call. = FALSE)
}

#' Sample size for average bioequivalence
#'
#' @param design `crossover` (two-period two-treatment) or `parallel`.
#' @param gmr Anticipated geometric mean ratio (test/reference).
#' @param cv_within Within-subject coefficient of variation on the raw
#'   (untransformed) scale. Required for `design = "crossover"`.
#' @param cv_between Between-subject coefficient of variation on the raw
#'   scale. Required for `design = "parallel"`.
#' @param lower Lower bioequivalence bound (ratio scale), usually `0.80`.
#' @param upper Upper bioequivalence bound (ratio scale), usually `1.25`.
#' @param alpha Type I error rate (one-sided; bioequivalence uses two
#'   one-sided tests, each at `alpha`).
#' @param power Target power.
#' @param allocation Allocation ratio `n_B / n_A`. Only used for
#'   `design = "parallel"`.
#' @param dropout Expected dropout proportion.
#' @param method `iterative_tost` (search for the smallest n achieving the
#'   target two-one-sided-test power, computed exactly via the noncentral t
#'   distribution; the preferred method) or `normal_approx` (closed-form
#'   approximation treating the variance as known).
#' @return A `sample_size` object. `n` is the raw (continuous) solution, as in
#'   the other `samplesize*()` functions. For the crossover design,
#'   `n_adjusted` and `n_total` are the recruitment target: the smallest even
#'   total (balanced sequences) that reaches the target power (checked with
#'   the exact TOST power for `iterative_tost`), inflated for dropout and
#'   rounded up again to an even number. For the parallel design,
#'   `n_adjusted`/`n_per_group` give the per-arm sizes and `n_total` their sum.
#' @details
#' On the log scale, with \eqn{\theta_0=\log(GMR)}, \eqn{\theta_L=\log(L)},
#' \eqn{\theta_U=\log(U)}, and \eqn{d_{BE}=\min(\theta_0-\theta_L,\
#' \theta_U-\theta_0)} (must be positive: the anticipated GMR must lie
#' strictly inside the bioequivalence bounds):
#'
#' `method = "iterative_tost"` searches for the smallest \eqn{n} such that
#' the exact two-one-sided-test power is at least the target. Let
#' \eqn{SE(n) = \sqrt{2\sigma_w^2/n}} for a crossover design (total \eqn{n},
#' \eqn{df=n-2}) or \eqn{SE(n_A) = \sigma_b\sqrt{(r+1)/(rn_A)}} for a
#' parallel design (per-arm \eqn{n_A}, \eqn{n_B=rn_A}, \eqn{df=n_A+n_B-2}).
#' The two one-sided test statistics \eqn{T_U=(\hat\theta_0-\theta_U)/SE}
#' and \eqn{T_L=(\hat\theta_0-\theta_L)/SE} each follow a *noncentral* t
#' distribution with the above \eqn{df} and noncentrality parameters
#' \eqn{\delta_U=(\theta_0-\theta_U)/SE}, \eqn{\delta_L=(\theta_0-\theta_L)/SE}
#' - not a location-shifted central t, and not a normal approximation, since
#' the estimated \eqn{SE} carries its own sampling variability. Power is
#' \deqn{Power(n) = P(T_U<-t_{1-\alpha,df}) + P(T_L>t_{1-\alpha,df}) - 1
#' = T_{df,\delta_U}(-t_{1-\alpha,df}) - T_{df,\delta_L}(t_{1-\alpha,df})}
#' where \eqn{T_{df,\delta}} is the noncentral t CDF (`stats::pt(...,
#' ncp = delta)`) - this is Phillips's (1990) formula, verified in package
#' tests to match `PowerTOST::power.TOST()` (the regulatory-standard,
#' Owen's-Q-based calculation) to numerical precision for balanced designs.
#' Coefficients of variation are converted to log-scale standard deviations
#' via \eqn{\sigma = \sqrt{\log(1+CV^2)}}.
#'
#' `method = "normal_approx"` uses the closed-form approximation
#' \deqn{n_{total} \approx \frac{2\sigma_w^2(z_{1-\beta}+z_{1-\alpha})^2}{d_{BE}^2}}
#' (crossover) or, per arm,
#' \deqn{n_A \approx \frac{(r+1)\sigma_b^2(z_{1-\beta}+z_{1-\alpha})^2}{rd_{BE}^2}}
#' (parallel), switching to \eqn{z_{1-\beta/2}} in place of \eqn{z_{1-\beta}}
#' when \eqn{GMR=1} exactly, matching the analogous special case in
#' [samplesizecontinuous()]. `iterative_tost` is preferred because this
#' approximation - which treats \eqn{\sigma} as known and uses normal
#' rather than noncentral-t quantiles - can under-size a study by a
#' material margin (verified up to ~20% at small n) away from \eqn{GMR=1}
#' or at small n.
#' @references
#' Julious SA. *Sample Sizes for Clinical Trials*. Chapman & Hall/CRC; 2010.
#'
#' Phillips KF. Power of the two one-sided tests procedure in
#' bioequivalence. *Journal of Pharmacokinetics and Biopharmaceutics*.
#' 1990;18(2):137-144.
#' @export
samplesizebioequivalence <- function(
  design = c("crossover", "parallel"),
  gmr = 1,
  cv_within = NULL,
  cv_between = NULL,
  lower = 0.80,
  upper = 1.25,
  alpha = 0.05,
  power = 0.90,
  allocation = 1,
  dropout = 0,
  method = c("iterative_tost", "normal_approx")
) {
  design <- match.arg(design)
  method <- match.arg(method)
  sample_size_validate_positive(gmr, "gmr")
  sample_size_validate_positive(lower, "lower")
  sample_size_validate_positive(upper, "upper")
  if (upper <= lower) {
    stop("`upper` must exceed `lower`.", call. = FALSE)
  }
  sample_size_validate_probability(alpha, "alpha")
  sample_size_validate_probability(power, "power")
  sample_size_validate_dropout(dropout)
  sample_size_validate_positive(allocation, "allocation")

  theta0 <- log(gmr)
  theta_l <- log(lower)
  theta_u <- log(upper)
  d_be <- min(theta0 - theta_l, theta_u - theta0)
  if (d_be <= 0) {
    stop("The anticipated GMR must lie strictly inside (lower, upper).", call. = FALSE)
  }
  z_alpha_one <- stats::qnorm(1 - alpha)
  z_power <- stats::qnorm(power)
  z_power_eq <- stats::qnorm(1 - (1 - power) / 2)
  at_center <- isTRUE(all.equal(theta0, 0))

  if (identical(design, "crossover")) {
    sample_size_validate_positive(cv_within, "cv_within")
    sigma <- sqrt(log(1 + cv_within^2))
    se_fn <- function(n) sqrt(2 * sigma^2 / n)
    df_fn <- function(n) n - 2
  } else {
    sample_size_validate_positive(cv_between, "cv_between")
    sigma <- sqrt(log(1 + cv_between^2))
    se_fn <- function(n) sigma * sqrt((allocation + 1) / (allocation * n))
    df_fn <- function(n) n * (1 + allocation) - 2
  }

  # Exact two-one-sided-test (TOST) power via the noncentral t distribution
  # (Phillips 1990, already cited below): with T_U = (theta0_hat -
  # theta_U)/SE_hat and T_L = (theta0_hat - theta_L)/SE_hat, each is
  # noncentral-t(df, ncp) under the true theta0 - not a location-shifted
  # central t, and not a normal approximation - and the joint TOST power is
  # bounded via P(T_U < -tcrit, T_L > tcrit) >= P(T_U < -tcrit) +
  # P(T_L > tcrit) - 1, which Phillips shows is close to exact in practice.
  # Verified against PowerTOST::power.TOST() (Owen's-Q-based, the
  # regulatory-standard exact calculation): matches to numerical precision
  # for balanced designs (see package tests); a prior version of this
  # function used stats::pnorm()/qnorm() (a normal approximation treating
  # sigma as known), which under-sized studies by up to ~20% at small n.
  tost_power <- function(n) {
    se <- se_fn(n)
    df <- df_fn(n)
    tcrit <- stats::qt(1 - alpha, df)
    delta_u <- (theta0 - theta_u) / se
    delta_l <- (theta0 - theta_l) / se
    stats::pt(-tcrit, df, ncp = delta_u) - stats::pt(tcrit, df, ncp = delta_l)
  }

  if (identical(method, "iterative_tost")) {
    n_raw <- stats::uniroot(function(n) tost_power(n) - power, c(4, 1e6), extendInt = "upX")$root
  } else if (identical(design, "crossover")) {
    n_raw <- if (at_center) {
      2 * sigma^2 * (z_power_eq + z_alpha_one)^2 / d_be^2
    } else {
      2 * sigma^2 * (z_power + z_alpha_one)^2 / d_be^2
    }
  } else {
    n_raw <- if (at_center) {
      (allocation + 1) * sigma^2 * (z_power_eq + z_alpha_one)^2 / (allocation * d_be^2)
    } else {
      (allocation + 1) * sigma^2 * (z_power + z_alpha_one)^2 / (allocation * d_be^2)
    }
  }

  method_label <- paste0(design, " bioequivalence (", method, ")")
  cv_note <- if (identical(design, "crossover")) {
    paste0("Within-subject CV = ", format_sample_value(cv_within), ".")
  } else {
    paste0("Between-subject CV = ", format_sample_value(cv_between), ".")
  }

  if (identical(design, "crossover")) {
    # A 2x2 crossover needs balanced sequences, so the total is even. With the
    # exact TOST power, step up from the rounded raw solution until the power
    # target is met (power is not perfectly continuous in n at small df).
    round_even <- function(n) as.integer(2 * ceiling(n / 2 - 1e-9))
    n_main <- round_even(n_raw)
    if (identical(method, "iterative_tost")) {
      while (tost_power(n_main) < power) n_main <- n_main + 2L
    }
    n_adj <- round_even(n_main / (1 - dropout))
    sample_size_result(
      endpoint = "bioequivalence", design = design, objective = "equivalence", method = method_label,
      n = n_raw, n_adjusted = n_adj, n_total = n_adj,
      assumptions = sample_size_assumptions(
        "A two-period two-treatment crossover with log-normal within-subject variation is assumed.",
        cv_note,
        paste0("Distance to the nearest bioequivalence boundary (log scale) = ", format_sample_value(d_be), ".")
      ),
      report = paste0(
        "Crossover bioequivalence sample size was calculated with gmr = ", format_sample_value(gmr),
        ", bounds = [", format_sample_value(lower), ", ", format_sample_value(upper), "]; required total n = ",
        format_sample_count(n_adj), "."
      ),
      plot_data = sample_size_plot_frame("Required total n", n_raw, n_adj)
    )
  } else {
    n_b_raw <- allocation * n_raw
    n_a_adj <- sample_size_apply_dropout(sample_size_round(n_raw), dropout)
    n_b_adj <- sample_size_apply_dropout(sample_size_round(n_b_raw), dropout)
    sample_size_result(
      endpoint = "bioequivalence", design = design, objective = "equivalence", method = method_label,
      n = n_raw, n_adjusted = c(A = n_a_adj, B = n_b_adj), n_per_group = c(A = n_a_adj, B = n_b_adj),
      n_total = n_a_adj + n_b_adj,
      assumptions = sample_size_assumptions(
        "Two independent groups with log-normal between-subject variation are assumed.",
        cv_note,
        paste0("Distance to the nearest bioequivalence boundary (log scale) = ", format_sample_value(d_be), ".")
      ),
      report = paste0(
        "Parallel bioequivalence sample size was calculated with gmr = ", format_sample_value(gmr),
        ", bounds = [", format_sample_value(lower), ", ", format_sample_value(upper), "]; required per-group sizes = A ",
        format_sample_count(n_a_adj), ", B ", format_sample_count(n_b_adj), "."
      ),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_raw, n_b_raw), c(n_a_adj, n_b_adj))
    )
  }
}

#' Sample size for ordinal endpoints
#'
#' @param design `parallel` only.
#' @param objective `superiority` only.
#' @param p_superiority Probability that a randomly selected subject in group A
#'   exceeds a subject in group B, with ties contributing half.
#' @param alpha Type I error rate.
#' @param power Target power.
#' @param dropout Expected dropout proportion.
#' @return A `sample_size` object.
#' @details
#' Noether's method for two independent ordinal (or continuous,
#' Mann-Whitney/Wilcoxon-type) groups, with
#' \eqn{P = P(A>B) + \tfrac12 P(A=B)} and equal allocation:
#' \deqn{n_{per\ group} = \frac{(z_{1-\alpha/2}+z_{1-\beta})^2}{6(P-0.5)^2}}
#' @references
#' Julious SA. *Sample Sizes for Clinical Trials*. Chapman & Hall/CRC; 2010.
#' @export
samplesizeordinal <- function(
  design = c("parallel"),
  objective = c("superiority"),
  p_superiority = NULL,
  alpha = 0.05,
  power = 0.90,
  dropout = 0
) {
  design <- match.arg(design)
  objective <- match.arg(objective)
  if (!identical(design, "parallel")) {
    stop("Only parallel ordinal sample size is implemented.", call. = FALSE)
  }
  if (!identical(objective, "superiority")) {
    stop("Ordinal sample size currently supports superiority only.", call. = FALSE)
  }
  sample_size_validate_probability(alpha, "alpha")
  sample_size_validate_probability(power, "power")
  sample_size_validate_dropout(dropout)
  z_power <- stats::qnorm(power)
  z_alpha_two <- stats::qnorm(1 - alpha / 2)

  sample_size_validate_probability(p_superiority, "p_superiority")
  if (p_superiority <= 0.5) {
    stop("p_superiority must exceed 0.5 for superiority planning.", call. = FALSE)
  }
  n_per_group_raw <- (z_alpha_two + z_power)^2 / (6 * (p_superiority - 0.5)^2)
  n_per_group <- sample_size_round(n_per_group_raw)
  n_adj <- sample_size_apply_dropout(n_per_group, dropout)
  sample_size_result(
    endpoint = "ordinal",
    design = "parallel",
    objective = objective,
    method = "Noether superiority approximation",
    n = n_per_group_raw,
    n_adjusted = c(A = n_adj, B = n_adj),
    n_per_group = c(A = n_adj, B = n_adj),
    n_total = 2 * n_adj,
    assumptions = sample_size_assumptions("Two independent ordinal groups are assumed."),
    report = paste0(
      "Ordinal superiority sample size was calculated with p_superiority = ",
      format_sample_value(p_superiority), "; required per-group size = ", format_sample_count(n_adj),
      ", total sample size = ", format_sample_count(2 * n_adj), "."
    ),
    plot_data = sample_size_plot_frame(c("A", "B"), c(n_per_group_raw, n_per_group_raw), c(n_adj, n_adj))
  )
}

#' Sample size for a correlation coefficient
#'
#' Number of subjects to detect a Pearson correlation `r` against `r0`,
#' using Fisher's z transformation.
#'
#' @param r Anticipated correlation.
#' @param r0 Correlation under the null hypothesis (default 0).
#' @param alpha Type I error rate.
#' @param power Target power.
#' @param alternative `"two.sided"` (default) or `"one.sided"`.
#' @param dropout Expected dropout proportion.
#' @return A `sample_size` object; `$n` is the unrounded n and
#'   `$n_adjusted` the rounded, dropout-inflated n.
#' @details
#' \deqn{n = \left(\frac{z_{1-\alpha/2} + z_{1-\beta}}{\operatorname{atanh} r - \operatorname{atanh} r_0}\right)^2 + 3}
#' (Rosner, Eq. 11.26), with \eqn{z_{1-\alpha}} for a one-sided test.
#' @references
#' Rosner B. *Fundamentals of Biostatistics* (8th ed.). Cengage; 2016.
#' Section 11.8.
#' @examples
#' samplesizecorrelation(r = 0.3)
#' @export
samplesizecorrelation <- function(r, r0 = 0, alpha = 0.05, power = 0.90, alternative = c("two.sided", "one.sided"),
                                  dropout = 0) {
  alternative <- match.arg(alternative)
  sample_size_validate_probability(alpha, "alpha"); sample_size_validate_probability(power, "power")
  sample_size_validate_dropout(dropout)
  if (abs(r) >= 1 || abs(r0) >= 1 || r == r0) stop("Need |r| < 1, |r0| < 1 and r != r0.", call. = FALSE)
  za <- stats::qnorm(1 - if (alternative == "two.sided") alpha / 2 else alpha)
  n_raw <- ((za + stats::qnorm(power)) / (atanh(r) - atanh(r0)))^2 + 3
  n_adj <- sample_size_apply_dropout(sample_size_round(n_raw), dropout)
  sample_size_result(
    endpoint = "correlation", design = "single sample", objective = "superiority",
    method = "Fisher z approximation", n = n_raw, n_adjusted = n_adj, n_total = n_adj,
    assumptions = sample_size_assumptions("Bivariate normal data are assumed.",
                                          paste0("Null correlation r0 = ", format_sample_value(r0), ".")),
    report = paste0("Correlation sample size was calculated with r = ", format_sample_value(r), ", r0 = ",
                    format_sample_value(r0), ", alpha = ", format_sample_value(alpha), " (", alternative, "), power = ",
                    format_sample_value(power), "; required sample size = ", format_sample_count(n_adj), "."),
    plot_data = sample_size_plot_frame("Required subjects", n_raw, n_adj))
}

#' Adjust a sample size for dropout
#'
#' @param n Evaluability sample size.
#' @param dropout Dropout proportion.
#' @return The dropout-adjusted sample size.
#' @export
samplesizeadjustdropout <- function(n, dropout = 0) {
  sample_size_validate_dropout(dropout)
  sample_size_round(n / (1 - dropout))
}

#' Adjust an individually randomized sample size for cluster randomization
#'
#' @param n_ind Individually randomized sample size (e.g. from
#'   [samplesizecontinuous()] or [samplesizebinary()]).
#' @param m Average cluster size.
#' @param rho Intracluster correlation.
#' @param cv_m Coefficient of variation of cluster sizes, for unequal
#'   cluster sizes. `NULL` (the default) assumes equal cluster sizes.
#' @return The cluster-adjusted sample size.
#' @details
#' Equal cluster size:
#' \deqn{DE = 1 + (m-1)\rho, \qquad n_{clustered} = \lceil n_{ind} \cdot DE \rceil}
#'
#' Unequal cluster size (approximate), with coefficient of variation
#' \eqn{CV_m}:
#' \deqn{DE \approx 1 + \left((1+CV_m^2)m - 1\right)\rho}
#' @references
#' Donner A, Klar N. *Design and Analysis of Cluster Randomization Trials in
#' Health Research*. Arnold; 2000.
#' @export
samplesizeclusteradjust <- function(n_ind, m, rho, cv_m = NULL) {
  sample_size_validate_positive(n_ind, "n_ind")
  sample_size_validate_positive(m, "m")
  sample_size_validate_probability(rho, "rho", allow_zero = TRUE)
  design_effect <- if (is.null(cv_m)) {
    1 + (m - 1) * rho
  } else {
    sample_size_validate_positive(cv_m, "cv_m")
    1 + ((1 + cv_m^2) * m - 1) * rho
  }
  sample_size_round(n_ind * design_effect)
}

#' Plot the cluster-adjusted sample size against the ICC
#'
#' Draws [samplesizeclusteradjust()] over a range of intracluster
#' correlations, one line per average cluster size, for a sensitivity figure
#' of a cluster randomized trial.
#'
#' @inheritParams samplesizeclusteradjust
#' @param m Average cluster size(s); one line per value.
#' @param rho Intracluster correlations on the x axis.
#' @return A `ggplot` object; its data (`rho`, `m`, `n`, `clusters`) hold the
#'   adjusted total size and the number of clusters, `ceiling(n / m)`.
#' @examples
#' plotclusteradjust(120, m = c(10, 20, 40), rho = seq(0, 0.1, by = 0.01))
#' @export
plotclusteradjust <- function(n_ind, m, rho = seq(0, 0.1, by = 0.01), cv_m = NULL) {
  d <- expand.grid(rho = rho, m = m)
  d$n <- mapply(function(r, k) samplesizeclusteradjust(n_ind, k, r, cv_m), d$rho, d$m); d$clusters <- ceiling(d$n / d$m)
  lev <- format(m); d$mf <- factor(format(d$m), lev)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$rho, y = .data$n, colour = .data$mf)) +
    ggplot2::geom_hline(yintercept = n_ind, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = tfpal(length(lev))) +
    ggplot2::labs(x = "Intracluster correlation", y = "Total sample size", colour = "Cluster size",
                  caption = paste0("Dashed: individually randomized size (", n_ind, "). Design effect ", if (is.null(cv_m))
                                     "1 + (m - 1) ICC." else paste0("1 + ((1 + CV^2) m - 1) ICC, CV = ", cv_m, "."))) +
    tftheme() + ggplot2::theme(legend.position = "top")
}

#' Sample size for confidence-interval precision
#'
#' @param endpoint Endpoint family: `continuous` or `binary`.
#' @param design `one_sample`, `two_sample`, or (binary only) `odds_ratio`.
#' @param width Desired confidence-interval half-width (margin of error), on
#'   the scale of the estimate (or of the log odds ratio for
#'   `design = "odds_ratio"`). This is NOT the full interval length: the
#'   planned interval is estimate +/- `width`, so its total length is
#'   `2 * width`. To plan for a total interval length `L`, pass
#'   `width = L / 2`.
#' @param sd Standard deviation (continuous endpoint).
#' @param p Anticipated proportion (binary, `one_sample`).
#' @param p1 Anticipated proportion in arm A (binary, `two_sample` or
#'   `odds_ratio`).
#' @param p2 Anticipated proportion in arm B (binary, `two_sample` or
#'   `odds_ratio`).
#' @param alpha Type I error rate (two-sided CI coverage is `1 - alpha`).
#' @param allocation Allocation ratio `n_B / n_A`. Used for
#'   `endpoint = "continuous"`, `design = "two_sample"`; binary
#'   `design = "two_sample"` (risk-difference half-width, unequal-allocation
#'   variance); and `design = "odds_ratio"`.
#' @param dropout Expected dropout proportion.
#' @param conservative For `endpoint = "binary"`, `design = "one_sample"`
#'   only: use the worst-case `p = 0.5` instead of the anticipated `p`. Only
#'   supported with `method = "wald"`.
#' @param method For `endpoint = "binary"`, `design = "one_sample"` only:
#'   the confidence-interval method used for planning. `"wald"` (default,
#'   backward-compatible) uses the closed-form normal-approximation formula.
#'   `"wilson"` and `"exact"` (Clopper-Pearson) instead run a deterministic
#'   integer search (see Details) for the minimum sample size whose
#'   confidence interval has a half-width no greater than `width`.
#' @param criterion For `endpoint = "binary"`, `design = "one_sample"`,
#'   `method %in% c("wilson", "exact")` only: `"expected"` (default)
#'   evaluates precision at the single anticipated event count
#'   `round(n * p)`. `"worst_case"` requires the precision target to hold
#'   over a prevalence-local band of plausible event counts (see Details).
#'   Recorded but does not alter the `"wald"` closed-form calculation.
#' @param min_expected_events For `endpoint = "binary"`, `design =
#'   "one_sample"` only: threshold below which the expected event count or
#'   expected non-event count triggers a rare-event diagnostic warning.
#' @param max_n For `endpoint = "binary"`, `design = "one_sample"`,
#'   `method %in% c("wilson", "exact")` only: upper bound for the minimum-n
#'   search. The search errors informatively if no `n <= max_n` satisfies
#'   the requested precision.
#' @return A `sample_size` object. For one-sample binary precision planning,
#'   `x$diagnostics` additionally holds rare-event and achieved-precision
#'   diagnostics (see Details).
#' @details
#' Unlike the other `sample_size_*()` functions, precision-based planning is
#' driven by a target confidence-interval half-width \eqn{w} rather than
#' power against an effect size.
#'
#' Continuous, one sample:
#' \deqn{n = \left(\frac{z_{1-\alpha/2}\sigma}{w}\right)^2}
#'
#' Continuous, two independent means:
#' \deqn{n_A = \frac{(r+1)\sigma^2z_{1-\alpha/2}^2}{rw^2}}
#'
#' Binary, one proportion, `method = "wald"` (normal approximation; the
#' original, backward-compatible formula):
#' \deqn{n = \frac{z_{1-\alpha/2}^2p(1-p)}{w^2}}
#' (or \eqn{n = z_{1-\alpha/2}^2/(4w^2)} for the conservative \eqn{p=0.5} case)
#'
#' Binary, one proportion, `method = "wilson"` or `"exact"`: there is no
#' closed form; \eqn{n} is found by a deterministic integer search. The
#' planning criterion is the maximum one-sided half-width
#' \eqn{\max(\hat p - lower,\ upper - \hat p) \le w}; because Wilson and
#' Clopper-Pearson intervals are asymmetric around \eqn{\hat p}, this is not
#' half of the total interval width. The interval is evaluated at:
#' \itemize{
#'   \item `method = "wilson"`, `criterion = "expected"`: the anticipated
#'     proportion itself, \eqn{\hat p = p}, with no rounding of the event
#'     count. The half-width then decreases steadily with \eqn{n}, and the
#'     result is the smallest \eqn{n} with
#'     \deqn{\frac{z}{1+z^2/n}\sqrt{\frac{p(1-p)}{n}+\frac{z^2}{4n^2}}
#'       + \left|\frac{p+z^2/(2n)}{1+z^2/n} - p\right| \le w,}
#'     which can be checked by hand.
#'   \item `method = "exact"`, `criterion = "expected"`: the anticipated
#'     integer count \eqn{x=\mathrm{round}(np)} (Clopper-Pearson limits are
#'     defined for integer counts).
#'   \item `criterion = "worst_case"`: every integer count in a
#'     prevalence-local band `qbinom(c(0.001, 0.999), n, p)`.
#' }
#' With integer counts the half-width is not monotone in \eqn{n}: it can meet
#' the target at some \eqn{n} and miss it again at a larger \eqn{n}. The
#' function then returns the smallest \eqn{n} from which the target is met
#' for every larger \eqn{n} (checked over a forward window of
#' \eqn{\max(50, 4/p)} values, and verified against a brute-force scan in
#' the package tests), so that recruiting a few more participants never
#' breaks the planned precision. A Clopper-Pearson interval guarantees
#' at-least-nominal coverage at the evaluated count(s); this is not a
#' guarantee on the achieved width at every conceivable outcome.
#'
#' A precision-based sample size answers "how wide will my confidence
#' interval be", not "how likely am I to observe any events at all"; for
#' rare-prevalence planning (e.g. newborn-screening or rare-adverse-event
#' studies) the expected event count, and the probability of observing zero
#' events, are reported as diagnostics precisely because a narrow interval
#' does not guarantee that any events will be observed.
#'
#' Binary, two independent proportions, risk-difference half-width, with
#' allocation ratio \eqn{r = n_B/n_A}:
#' \deqn{n_A = \frac{z_{1-\alpha/2}^2}{w^2}\left[p_A(1-p_A) + p_B(1-p_B)/r\right], \qquad n_B = r n_A}
#' (reduces to the earlier equal-allocation-only formula when \eqn{r=1}).
#'
#' Binary, log-odds-ratio half-width:
#' \deqn{n_A = \frac{z_{1-\alpha/2}^2}{w^2}
#' \left[\frac1{p_A}+\frac1{1-p_A}+\frac1{rp_B}+\frac1{r(1-p_B)}\right]}
#' @references
#' Julious SA. *Sample Sizes for Clinical Trials*. Chapman & Hall/CRC; 2010.
#'
#' Wilson EB. Probable inference, the law of succession, and statistical
#' inference. *Journal of the American Statistical Association*.
#' 1927;22(158):209-212.
#'
#' Clopper CJ, Pearson ES. The use of confidence or fiducial limits
#' illustrated in the case of the binomial. *Biometrika*. 1934;26(4):404-413.
#' @examples
#' # 1) Common prevalence, Wald (backward-compatible default)
#' samplesizeprecision(endpoint = "binary", design = "one_sample", width = 0.05, p = 0.30)
#'
#' # 2) Rare prevalence, Wilson score, iterative search
#' samplesizeprecision(
#'   endpoint = "binary", design = "one_sample", width = 0.02,
#'   p = 0.01, method = "wilson"
#' )
#'
#' # 3) Rare prevalence, exact Clopper-Pearson, worst-case criterion
#' samplesizeprecision(
#'   endpoint = "binary", design = "one_sample", width = 0.005,
#'   p = 0.001, method = "exact", criterion = "worst_case"
#' )
#'
#' # 4) Dropout inflation: complete_case_n vs. adjusted_n
#' x <- samplesizeprecision(
#'   endpoint = "binary", design = "one_sample", width = 0.02,
#'   p = 0.05, method = "wilson", dropout = 0.10
#' )
#' x$diagnostics$complete_case_n
#' x$diagnostics$adjusted_n
#'
#' # 5) Two independent proportions, unequal allocation
#' samplesizeprecision(
#'   endpoint = "binary", design = "two_sample", width = 0.08,
#'   p1 = 0.4, p2 = 0.3, allocation = 2
#' )
#' @export
samplesizeprecision <- function(
  endpoint = c("continuous", "binary"),
  design = c("one_sample", "two_sample", "odds_ratio"),
  width,
  sd = NULL,
  p = NULL,
  p1 = NULL,
  p2 = NULL,
  alpha = 0.05,
  allocation = 1,
  dropout = 0,
  conservative = FALSE,
  method = c("wald", "wilson", "exact"),
  criterion = c("expected", "worst_case"),
  min_expected_events = 5,
  max_n = 1e7
) {
  endpoint <- match.arg(endpoint)
  design <- match.arg(design)
  sample_size_validate_probability(alpha, "alpha")
  sample_size_validate_dropout(dropout)
  sample_size_validate_positive(allocation, "allocation")
  sample_size_validate_positive(width, "width")
  if (identical(endpoint, "continuous") && identical(design, "odds_ratio")) {
    stop("`design = \"odds_ratio\"` is only available for `endpoint = \"binary\"`.", call. = FALSE)
  }
  z_alpha <- stats::qnorm(1 - alpha / 2)

  if (identical(endpoint, "continuous")) {
    sample_size_validate_positive(sd, "sd")
    if (identical(design, "one_sample")) {
      n_raw <- (z_alpha * sd / width)^2
      n_adj <- sample_size_apply_dropout(sample_size_round(n_raw), dropout)
      return(sample_size_result(
        endpoint = endpoint, design = design, objective = "precision",
        method = "one-sample CI half-width",
        n = n_raw, n_adjusted = n_adj, n_total = n_adj,
        assumptions = sample_size_assumptions(paste0("Target CI half-width = ", format_sample_value(width), ".")),
        report = paste0("One-sample precision sample size was calculated with sd = ", format_sample_value(sd), ", width = ", format_sample_value(width), "; required n = ", format_sample_count(n_adj), "."),
        plot_data = sample_size_plot_frame("Required n", n_raw, n_adj)
      ))
    }
    n_a_raw <- ((allocation + 1) * sd^2 * z_alpha^2) / (allocation * width^2)
    n_b_raw <- allocation * n_a_raw
    n_a_adj <- sample_size_apply_dropout(sample_size_round(n_a_raw), dropout)
    n_b_adj <- sample_size_apply_dropout(sample_size_round(n_b_raw), dropout)
    return(sample_size_result(
      endpoint = endpoint, design = design, objective = "precision",
      method = "two-sample CI half-width",
      n = n_a_raw, n_adjusted = c(A = n_a_adj, B = n_b_adj), n_per_group = c(A = n_a_adj, B = n_b_adj),
      assumptions = sample_size_assumptions(paste0("Target CI half-width for the mean difference = ", format_sample_value(width), ".")),
      report = paste0("Two-sample precision sample size was calculated with sd = ", format_sample_value(sd), ", width = ", format_sample_value(width), "; required per-group sizes = A ", format_sample_count(n_a_adj), ", B ", format_sample_count(n_b_adj), "."),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_a_raw, n_b_raw), c(n_a_adj, n_b_adj))
    ))
  }

  if (identical(design, "one_sample")) {
    return(sample_size_precision_binary_one_sample(
      width = width, p = p, alpha = alpha, dropout = dropout, conservative = conservative,
      method = method, criterion = criterion, min_expected_events = min_expected_events, max_n = max_n,
      z_alpha = z_alpha
    ))
  }

  sample_size_validate_probability(p1, "p1")
  sample_size_validate_probability(p2, "p2")
  if (identical(design, "two_sample")) {
    n_a_raw <- (z_alpha^2 / width^2) * (p1 * (1 - p1) + p2 * (1 - p2) / allocation)
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = endpoint, design = design, objective = "precision",
      method = "two-proportion CI half-width (unequal-allocation risk-difference variance)",
      n = n_a_raw, n_adjusted = c(A = n_a_adj, B = n_b_adj), n_per_group = c(A = n_a_adj, B = n_b_adj), n_total = n_a_adj + n_b_adj,
      assumptions = sample_size_assumptions(
        paste0("Target CI half-width for the risk difference = ", format_sample_value(width), "."),
        paste0("Allocation ratio n_B / n_A = ", format_sample_value(allocation), "; Var(p1_hat - p2_hat) = p1(1-p1)/n_A + p2(1-p2)/n_B is used directly (this generalizes the equal-allocation formula, which is the r = 1 special case).")
      ),
      report = paste0("Two-proportion precision sample size was calculated with p1 = ", format_sample_value(p1), ", p2 = ", format_sample_value(p2), ", width = ", format_sample_value(width), ", allocation ratio = ", format_sample_value(allocation), "; required per-group sizes = A ", format_sample_count(n_a_adj), ", B ", format_sample_count(n_b_adj), "."),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_a_raw, n_b_raw), c(n_a_adj, n_b_adj))
    ))
  }

  n_a_raw <- (z_alpha^2 / width^2) * (1 / p1 + 1 / (1 - p1) + 1 / (allocation * p2) + 1 / (allocation * (1 - p2)))
  n_b_raw <- allocation * n_a_raw
  n_a_adj <- sample_size_apply_dropout(sample_size_round(n_a_raw), dropout)
  n_b_adj <- sample_size_apply_dropout(sample_size_round(n_b_raw), dropout)
  sample_size_result(
    endpoint = endpoint, design = design, objective = "precision",
    method = "log-odds-ratio CI half-width",
    n = n_a_raw, n_adjusted = c(A = n_a_adj, B = n_b_adj), n_per_group = c(A = n_a_adj, B = n_b_adj),
    assumptions = sample_size_assumptions(paste0("Target CI half-width on the log-odds-ratio scale = ", format_sample_value(width), ".")),
    report = paste0("Log-odds-ratio precision sample size was calculated with p1 = ", format_sample_value(p1), ", p2 = ", format_sample_value(p2), ", width = ", format_sample_value(width), "; required per-group sizes = A ", format_sample_count(n_a_adj), ", B ", format_sample_count(n_b_adj), "."),
    plot_data = sample_size_plot_frame(c("A", "B"), c(n_a_raw, n_b_raw), c(n_a_adj, n_b_adj))
  )
}

# --- One-proportion precision planning (endpoint = "binary", design = "one_sample") ---
#
# Dispatches across method = wald/wilson/exact and criterion =
# expected/worst_case, builds rare-event and achieved-precision diagnostics,
# and assembles the sample_size_result(). Kept out of samplesizeprecision()
# itself to keep that dispatcher flat.
sample_size_precision_binary_one_sample <- function(width, p, alpha, dropout, conservative, method, criterion, min_expected_events, max_n, z_alpha) {
  if (width >= 1) {
    stop("`width` must be in (0, 1) for binary precision planning.", call. = FALSE)
  }
  method <- match.arg(method, c("wald", "wilson", "exact"))
  criterion <- match.arg(criterion, c("expected", "worst_case"))
  sample_size_validate_nonnegative(min_expected_events, "min_expected_events")
  max_n <- sample_size_validate_max_n(max_n)

  if (conservative && !identical(method, "wald")) {
    stop("`conservative = TRUE` is only supported for `method = \"wald\"`; Wilson and exact planning require an anticipated `p` (there is no worst-case p for an asymmetric search).", call. = FALSE)
  }

  if (conservative) {
    p_effective <- 0.5
    n_raw <- z_alpha^2 / (4 * width^2)
    complete_case_n <- sample_size_round(n_raw)
  } else {
    sample_size_validate_probability(p, "p")
    p_effective <- p
    if (identical(method, "wald")) {
      n_raw <- z_alpha^2 * p * (1 - p) / width^2
      complete_case_n <- sample_size_round(n_raw)
    } else {
      complete_case_n <- sample_size_binomial_required_n(p, width, alpha, method, criterion, max_n)
      n_raw <- complete_case_n
    }
  }
  n_adj <- sample_size_apply_dropout(complete_case_n, dropout)

  diagnostics <- sample_size_binary_precision_diagnostics(
    p = p_effective, complete_case_n = complete_case_n, n_adj = n_adj, width = width,
    alpha = alpha, method = method, criterion = criterion, min_expected_events = min_expected_events,
    dropout = dropout
  )

  ci_display_name <- switch(method,
    wald = "normal-approximation (Wald)",
    wilson = "Wilson score",
    exact = "Clopper-Pearson exact"
  )
  method_label <- paste0(
    "one-proportion CI half-width (", method, if (conservative) ", conservative p = 0.5" else "", ", ", criterion, " criterion)"
  )
  p_label <- if (conservative) "conservative p = 0.5" else paste0("p = ", format_sample_value(p_effective))
  conf_pct <- formatC(100 * (1 - alpha), format = "f", digits = 0)
  dropout_pct <- formatC(100 * dropout, format = "f", digits = 0)

  report <- paste0(
    "One-proportion precision planning used a ", conf_pct, "% ", ci_display_name, " confidence interval, ",
    if (conservative) "a conservative p = 0.5" else paste0("anticipated prevalence p = ", format_sample_value(p_effective)),
    ", and requested maximum half-width ", format_sample_value(width), ". ",
    "The minimum complete-case sample size was ", format_sample_count(complete_case_n), ". ",
    if (dropout > 0) {
      paste0("After ", dropout_pct, "% dropout inflation, the recruitment target was ", format_sample_count(n_adj), ". ")
    } else {
      paste0("No dropout inflation was applied; the recruitment target was ", format_sample_count(n_adj), ". ")
    },
    "The expected number of events is ", format_sample_value(diagnostics$expected_events),
    " and the achieved maximum half-width is ", format_sample_value(diagnostics$achieved_maximum_half_width), "."
  )

  sample_size_result(
    endpoint = "binary", design = "one_sample", objective = "precision",
    method = method_label,
    n = n_raw, n_adjusted = n_adj, n_total = n_adj,
    assumptions = sample_size_assumptions(
      paste0("Confidence interval method: ", ci_display_name, "."),
      if (identical(criterion, "expected")) {
        if (identical(method, "wilson")) "The interval is evaluated at the anticipated proportion p itself (no rounding of the event count)." else "Event counts are evaluated at the single anticipated event count round(n * p)."
      } else {
        "Event counts are evaluated over a prevalence-local band of plausible counts (binomial 0.001/0.999 quantiles around n * p), not the anticipated count alone."
      },
      "\"Half-width\" is the maximum one-sided distance from the reference proportion to either confidence limit; for asymmetric (Wilson/exact) intervals this is not the same as half of the total interval width.",
      if (identical(method, "wald")) {
        "Wald planning uses a closed-form normal approximation and can be unreliable for small or extreme proportions."
      } else {
        "Wilson/exact planning uses a deterministic integer numerical search for the minimum sample size (no closed-form formula exists)."
      },
      if (identical(method, "exact")) {
        "Clopper-Pearson exact planning guarantees at-least-nominal coverage at the evaluated event count(s); this is not a guarantee on the achieved width at every conceivable outcome."
      } else {
        NULL
      },
      "Precision-based planning targets confidence-interval width, not the probability of observing any events; a narrow planned interval does not by itself ensure that cases will be observed (see the rare-event diagnostics)."
    ),
    report = report,
    plot_data = sample_size_plot_frame("Required n", complete_case_n, n_adj),
    diagnostics = diagnostics
  )
}

# Rare-event and achieved-precision diagnostics for one-proportion precision
# planning. Computes the achieved interval at the anticipated event count
# round(complete_case_n * p) using the same CI method selected for planning;
# for criterion = "worst_case" (non-Wald), achieved_maximum_half_width instead
# reports the actual worst-case value used by the search. Any warnings are
# collected and issued exactly once here (never inside the search helpers).
#
# For method = "wald", the reference proportion is p itself, not a rounded
# event count: the closed-form n_raw = z^2*p*(1-p)/w^2 already guarantees
# z*sqrt(p(1-p)/n) <= w at continuous p, and complete_case_n = ceiling(n_raw)
# preserves that; using round(n*p)/n instead would reintroduce integer
# rounding noise the closed form was never meant to have.
sample_size_binary_precision_diagnostics <- function(p, complete_case_n, n_adj, width, alpha, method, criterion, min_expected_events, dropout) {
  if (identical(method, "wald") || sample_size_binomial_continuous(method, criterion)) {
    x_expected <- p * complete_case_n
    p_ref <- p
  } else {
    x_expected <- round(complete_case_n * p)
    p_ref <- x_expected / complete_case_n
  }
  ci <- sample_size_binomial_half_width(x_expected, complete_case_n, alpha, method, p_ref)
  achieved_maximum_half_width <- if (!identical(method, "wald")) {
    sample_size_binomial_max_half_width(complete_case_n, p, alpha, method, criterion)
  } else {
    ci$maximum_half_width
  }

  expected_events <- n_adj * p
  expected_non_events <- n_adj * (1 - p)

  messages <- character(0)
  if (expected_events < min_expected_events || expected_non_events < min_expected_events) {
    messages <- c(messages, paste0(
      "Expected event or non-event count is below min_expected_events = ", min_expected_events,
      "; rare-event asymptotics may be unreliable regardless of the confidence-interval method."
    ))
  }
  if (identical(method, "wald")) {
    messages <- c(messages, "Wald normal-approximation intervals can be unreliable for small or extreme proportions; consider method = \"wilson\" or method = \"exact\".")
  }
  approximation_warning <- if (length(messages)) paste(messages, collapse = " ") else NULL
  if (!is.null(approximation_warning)) warning(approximation_warning, call. = FALSE)

  list(
    anticipated_prevalence = p,
    expected_events_raw = complete_case_n * p,
    expected_events = expected_events,
    expected_non_events = expected_non_events,
    probability_zero_events = exp(n_adj * log1p(-p)),
    probability_at_least_one_event = -expm1(n_adj * log1p(-p)),
    requested_half_width = width,
    achieved_lower = ci$lower,
    achieved_upper = ci$upper,
    achieved_lower_half_width = ci$lower_half_width,
    achieved_upper_half_width = ci$upper_half_width,
    achieved_maximum_half_width = achieved_maximum_half_width,
    achieved_total_width = ci$total_width,
    confidence_level = 1 - alpha,
    ci_method = method,
    precision_criterion = criterion,
    complete_case_n = complete_case_n,
    adjusted_n = n_adj,
    dropout = dropout,
    approximation_warning = approximation_warning
  )
}

# --- Binomial confidence-interval helpers (base R only) ---
#
# Shared by sample_size_precision_binary_one_sample() for both diagnostics on
# a fixed n and the minimum-n search. `method` is one of "wald", "wilson",
# "exact"; intervals are truncated to [0, 1] after computation.

sample_size_binomial_ci <- function(x, n, alpha, method) {
  z <- stats::qnorm(1 - alpha / 2)
  p_hat <- x / n
  if (identical(method, "wald")) {
    se <- sqrt(p_hat * (1 - p_hat) / n)
    lower <- p_hat - z * se
    upper <- p_hat + z * se
  } else if (identical(method, "wilson")) {
    denom <- 1 + z^2 / n
    center <- (p_hat + z^2 / (2 * n)) / denom
    half <- (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
    lower <- center - half
    upper <- center + half
  } else {
    lower <- if (x == 0) 0 else stats::qbeta(alpha / 2, x, n - x + 1)
    upper <- if (x == n) 1 else stats::qbeta(1 - alpha / 2, x + 1, n - x)
  }
  list(lower = max(0, min(1, lower)), upper = max(0, min(1, upper)))
}

# Half-width relative to a reference proportion (not necessarily p_hat = x/n):
# for asymmetric intervals, the one-sided distances to each limit differ, and
# the planning criterion is their maximum, never half of the total width.
sample_size_binomial_half_width <- function(x, n, alpha, method, p_reference) {
  ci <- sample_size_binomial_ci(x, n, alpha, method)
  lower_half_width <- p_reference - ci$lower
  upper_half_width <- ci$upper - p_reference
  list(
    lower = ci$lower,
    upper = ci$upper,
    lower_half_width = lower_half_width,
    upper_half_width = upper_half_width,
    maximum_half_width = max(lower_half_width, upper_half_width),
    total_width = ci$upper - ci$lower
  )
}

# Candidate event counts for a given planning criterion at sample size n.
# "expected": the single anticipated count round(n * p).
# "worst_case": a prevalence-local band of counts (binomial 0.001/0.999
# quantiles) - far cheaper and more relevant to the planned p than a global
# search over x = 0..n, which would mostly probe implausible outcomes. When
# the band exceeds 10,000 integers (only for very large n), a deterministic
# grid (both endpoints, floor/ceiling/round of n*p, and an evenly spaced
# 2,000-point subset) approximates it instead of enumerating every integer.
sample_size_binomial_candidate_counts <- function(n, p, criterion, method = "exact") {
  if (identical(criterion, "expected")) {
    if (sample_size_binomial_continuous(method, criterion)) return(n * p)
    return(round(n * p))
  }
  q <- stats::qbinom(c(0.001, 0.999), size = n, prob = p)
  lo <- max(0, q[1])
  hi <- min(n, q[2])
  if (hi - lo + 1 > 10000) {
    grid <- c(lo, hi, floor(n * p), ceiling(n * p), round(n * p), as.integer(round(seq(lo, hi, length.out = 2000))))
    grid <- grid[grid >= 0 & grid <= n]
    return(sort(unique(grid)))
  }
  lo:hi
}

# Wilson (and Wald) intervals are defined for a non-integer count, so the
# "expected" criterion evaluates them at the anticipated proportion p itself
# (x = n * p). This removes the rounding noise of round(n * p) and makes the
# criterion monotone in n. Clopper-Pearson limits need integer counts.
sample_size_binomial_continuous <- function(method, criterion) {
  identical(criterion, "expected") && method %in% c("wilson", "wald")
}

sample_size_binomial_max_half_width <- function(n, p, alpha, method, criterion) {
  counts <- sample_size_binomial_candidate_counts(n, p, criterion, method)
  hws <- vapply(counts, function(x) {
    sample_size_binomial_half_width(x, n, alpha, method, p_reference = x / n)$maximum_half_width
  }, numeric(1))
  max(hws)
}

# Deterministic search for Wilson/exact planning (no simulation, no random
# numbers). Returns the smallest n from which the precision target holds for
# every larger n ("stable" n), not merely the first n that happens to pass.
#
# Why: with integer event counts (method = "exact", or criterion =
# "worst_case"), max_half_width(n) is not monotone in n. It can pass at n and
# fail again at n + k when round(n * p) or the qbinom() band jumps, so the
# first passing n would be a fragile plan: recruiting a few more
# participants would break the target. For Wilson with criterion =
# "expected" the count is continuous (n * p) and the criterion is monotone;
# the same algorithm then returns the plain minimum.
#
# Steps: (1) doubling from the Wald closed form until some n passes, then
# bisection to a passing boundary n; (2) forward check over a window of
# max(50, 4 / p) values above n (spans several round(n * p) jumps); on any
# failure, move n just past the last failure and check again, until a full
# window is clean; (3) walk backward while n - 1 still passes, so the result
# is the smallest n of the passing run. Verified against a brute-force scan
# in the package tests.
sample_size_binomial_search <- function(p, width, alpha, method, criterion, max_n) {
  passes <- function(n) sample_size_binomial_max_half_width(n, p, alpha, method, criterion) <= width
  no_solution <- function() stop(
    "No sample size up to max_n = ", max_n, " achieves the requested precision ",
    "(method = \"", method, "\", criterion = \"", criterion, "\"). ",
    "Increase `max_n` or relax `width`.",
    call. = FALSE
  )

  z <- stats::qnorm(1 - alpha / 2)
  hi <- min(max(2, ceiling(z^2 * p * (1 - p) / width^2)), max_n)
  while (!passes(hi)) {
    if (hi >= max_n) no_solution()
    hi <- min(hi * 2, max_n)
  }
  lo <- 1L
  while (hi - lo > 1) {
    mid <- lo + (hi - lo) %/% 2L
    if (passes(mid)) hi <- mid else lo <- mid
  }

  window <- max(50L, as.integer(ceiling(4 / p)))
  repeat {
    upper <- min(max_n, hi + window)
    if (upper <= hi) break
    candidates <- (hi + 1L):upper
    ok <- vapply(candidates, passes, logical(1))
    if (all(ok)) break
    hi <- max(candidates[!ok]) + 1L
    if (hi > max_n) no_solution()
    while (!passes(hi)) {
      hi <- hi + 1L
      if (hi > max_n) no_solution()
    }
  }

  while (hi > 2L && passes(hi - 1L)) hi <- hi - 1L
  as.integer(hi)
}

sample_size_binomial_required_n <- function(p, width, alpha, method, criterion, max_n) {
  sample_size_binomial_search(p, width, alpha, method, criterion, max_n)
}

#' Plot a sample size object
#'
#' @param x A sample-size object returned by `samplesize()`.
#' @param type Plot style: `both`, `curve`, or `summary`.
#' @param vary Optional named list with one planning input of the design and
#'   the values to try, for example `list(delta = c(2, 3, 4))`, `list(p1 =
#'   c(0.25, 0.30))` or `list(hr = c(0.6, 0.7, 0.8))`. The design is rerun for
#'   each value and one power curve per value is drawn against the total size
#'   before dropout, with the size needed for the planned power marked. Power is
#'   exact (noncentral t) for `method = "t"` and follows the design's own
#'   formula otherwise. Default `NULL` draws the single planned curve. For
#'   survival designs, varying `hr` recomputes the experimental-arm survival as
#'   `survival_a^hr` (proportional hazards) when `survival_a` is given.
#' @param what With `vary`, the x axis: `"n"` (total patients, default) or
#'   `"events"` (total events, survival designs only).
#' @param varylabels With `vary`, legend text for the values: a character
#'   vector in the order of the values, or named by the values (for example
#'   `c("0.7" = "HR 0.7")`). Default `NULL` prints the values to 3 significant
#'   digits.
#' @param ... Unused.
#' @details The curve shows power against sample size. Designs without an
#'   exact curve use the normal approximation scaled to the computed size:
#'   since \eqn{n \propto (z_{1-lpha/2} + z_{1-eta})^2}, power at \eqn{n}
#'   is \eqn{\Phi(\sqrt{n / n_0}(z_{1-lpha/2} + z_{1-eta}) - z_{1-lpha/2})},
#'   with \eqn{n_0} the unrounded total. Designs without a power target
#'   (precision) show the summary bars only. Under `method = "t"` the curve
#'   is the exact noncentral-t power.
#' @return A `ggplot2` object, or a combined patchwork object when `type =
#'   "both"`.
#' @export
plot.sample_size <- function(x, type = c("both", "curve", "summary"), vary = NULL, what = c("n", "events"),
                             varylabels = NULL, ...) {
  type <- match.arg(type)
  if (!is.null(vary)) return(ssfamily(x, vary, match.arg(what), varylabels))
  summary_plot <- NULL
  curve_plot <- NULL
  pd <- x$plot_data; sub <- paste0(x$objective, ", ", x$method)

  if (type != "curve" && !is.null(pd) && nrow(pd)) {
    # Raw sizes are fractional; label them to one decimal so the bar and its label agree.
    pd$display <- ifelse(pd$type == "raw", formatC(pd$value, format = "f", digits = 1), pd$display)
    pd$type <- factor(pd$type, c("raw", "adjusted"), c("Raw", "Planned"))
    summary_plot <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$label, y = .data$value, fill = .data$type)) +
      ggplot2::geom_col(width = 0.7, colour = "white", position = ggplot2::position_dodge(width = 0.75)) +
      ggplot2::geom_text(ggplot2::aes(label = .data$display), vjust = -0.35, size = 3.3,
                         position = ggplot2::position_dodge(width = 0.75)) +
      ggplot2::scale_fill_manual(values = c("grey75", TFCOL[["a1"]])) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0, 0.12))) +
      ggplot2::labs(subtitle = sub, x = NULL, y = "Sample size",
                    fill = NULL) +
      tftheme() + ggplot2::theme(legend.position = "top", panel.grid.major.x = ggplot2::element_blank())
  }

  cd <- x$curve_data
  pd <- x$plot_data
  if (!is.null(x$powerfun)) {
    tn <- sum(pd$value[pd$type == "adjusted"]); ng <- sample_size_curve_grid(tn, length.out = 60)
    cd <- sample_size_curve_bundle(tn, x$power, ng, vapply(ng, x$powerfun, 0), sslab(pd), "noncentral t")
  }
  if ((is.null(cd) || !nrow(cd)) && length(x$power) == 1 && length(x$alpha) == 1 && !is.null(pd) && nrow(pd)) {
    # Normal approximation: n is proportional to (z_alpha + z_beta)^2, so power at n follows from the computed n.
    raw <- sum(pd$value[pd$type == "raw"]); adj <- sum(pd$value[pd$type == "adjusted"])
    za <- stats::qnorm(1 - x$alpha / 2); eq <- identical(x$objective, "equivalence")
    if (eq) za <- stats::qnorm(1 - x$alpha)
    zb <- stats::qnorm(if (eq) (1 + x$power) / 2 else x$power)
    ng <- sample_size_curve_grid(adj, length.out = 60)
    pw <- stats::pnorm(sqrt(ng / raw) * (za + zb) - za); if (eq) pw <- pmax(2 * pw - 1, 0)
    cd <- sample_size_curve_bundle(adj, x$power, ng, pw, sslab(pd), "normal approximation")
  }
  if (type != "summary" && !is.null(cd) && is.data.frame(cd) && nrow(cd)) {
    tn <- cd$target_n[1]; pt <- cd$power_target[1]
    curve_plot <- ggplot2::ggplot(cd, ggplot2::aes(x = .data$n, y = .data$power)) +
      ggplot2::geom_hline(yintercept = pt, linetype = "dashed", colour = "grey55") +
      ggplot2::geom_vline(xintercept = tn, linetype = "dashed", colour = "grey55") +
      ggplot2::geom_line(colour = TFCOL[["a1"]], linewidth = 0.9) +
      ggplot2::annotate("point", x = tn, y = pt, colour = TFCOL[["a2"]], size = 2.5) +
      ggplot2::annotate("text", x = tn, y = 0.05, hjust = -0.08, vjust = 0, size = 3,
                        label = paste0("planned n = ", format_sample_count(tn))) +
      ggplot2::scale_y_continuous(labels = function(b) paste0(100 * b, "%"), limits = c(0, 1)) +
      ggplot2::labs(subtitle = sub, caption = switch(cd$power_label[1],
                      "normal approximation" = "Power curve: normal approximation scaled to the computed sample size.",
                      "noncentral t" = "Power curve: exact noncentral t."),
                    x = cd$n_label[1], y = "Power") +
      tftheme()
  }

  if (type == "summary") return(summary_plot)
  if (type == "curve") return(curve_plot %||% summary_plot)
  if (!is.null(curve_plot) && !is.null(summary_plot)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(patchwork::wrap_plots(curve_plot, summary_plot + ggplot2::labs(subtitle = NULL), ncol = 2,
                                   widths = c(3, 2)))
    }
    return(curve_plot)
  }
  curve_plot %||% summary_plot
}

# x-axis label of the power curve: total size for two arms, else the planned unit ("Pairs", "Total sample size").
sslab <- function(pd) {
  lb <- if (sum(pd$type == "raw") > 1) "Total sample size" else sub("^Required ", "", pd$label[1])
  paste0(toupper(substr(lb, 1, 1)), substring(lb, 2))
}

# One power curve per value of a planning input: rerun the design at that value over a grid of target powers
# (each run gives a point (n, power) on the design's own curve), or use its exact power function.
ssfamily <- function(x, vary, what = "n", varylabels = NULL) {
  if (is.null(x$fun) || is.null(x$power) || length(vary) != 1 || is.null(names(vary)))
    stop("`vary` needs a named list with one planning input, on a design with a power target.", call. = FALSE)
  ev <- what == "events"
  if (ev && is.null(x$required_events)) stop("`what = \"events\"` needs a survival design.", call. = FALSE)
  nm <- names(vary); a <- x$inputs; if ("dropout" %in% names(a)) a$dropout <- 0
  # Under proportional hazards the experimental-arm survival moves with the hazard ratio.
  ph <- nm == "hr" && !is.null(a$survival_a)
  run <- function(v, p) { a[[nm]] <- v; a$power <- p; if (ph) a$survival_b <- a$survival_a^v; do.call(x$fun, a) }
  tot <- function(r) if (ev) r$required_events else sum(r$plot_data$value[r$plot_data$type == "raw"])
  rs <- lapply(vary[[1]], run, p = x$power); tn <- vapply(rs, tot, 0); xmax <- 1.35 * max(tn)
  val <- as.character(signif(vary[[1]], 3)); ps <- seq(0.1, 0.99, by = 0.01)
  lev <- if (is.null(varylabels)) val else if (is.null(names(varylabels))) as.character(varylabels) else
    ifelse(val %in% names(varylabels), varylabels[val], val)
  cv <- do.call(rbind, lapply(seq_along(rs), function(i) {
    d <- if (!is.null(rs[[i]]$powerfun)) {
      ng <- seq(2, xmax, length.out = 80); data.frame(n = ng, power = vapply(ng, rs[[i]]$powerfun, 0))
    } else data.frame(n = vapply(ps, function(p) tot(run(vary[[1]][i], p)), 0), power = ps)
    d$v <- lev[i]; d[d$n <= xmax, ]
  }))
  cv$v <- factor(cv$v, lev); tg <- data.frame(n = tn, power = x$power, v = factor(lev, lev))
  ex <- !is.null(rs[[1]]$powerfun)
  ggplot2::ggplot(cv, ggplot2::aes(x = .data$n, y = .data$power, colour = .data$v)) +
    ggplot2::geom_hline(yintercept = x$power, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(data = tg, size = 2.3) +
    ggplot2::geom_text(ggplot2::aes(label = format(ceiling(.data$n), big.mark = ",")), data = tg, vjust = -0.9, size = 3,
                       show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = tfpal(length(lev))) +
    ggplot2::scale_y_continuous(labels = function(b) paste0(100 * b, "%"), limits = c(0, 1)) +
    ggplot2::scale_x_continuous(limits = c(0, xmax)) +
    ggplot2::labs(subtitle = paste0(x$objective, ", ", x$method),
                  x = if (ev) "Total events" else paste0(sslab(x$plot_data), ", before dropout"),
                  y = "Power", colour = nm, caption = paste0(
                    if (ex) "Exact noncentral t power. " else "Power from the design formula. ",
                    "Points: size for ", 100 * x$power, "% power.")) +
    tftheme() + ggplot2::theme(legend.position = "top")
}

#' Print a sample size object
#'
#' @param x A sample-size object.
#' @param ... Unused.
#' @return The object invisibly.
#' @export
print.sample_size <- function(x, ...) {
  tf_with_cli_colors({
    tf_title("Sample size planning")
    tf_field("Endpoint", x$endpoint)
    tf_field("Design", x$design)
    tf_field("Objective", x$objective)
    tf_field("Method", x$method)
    tf_blank()
    if (!is.null(x$assumptions)) {
      tf_section("Assumptions")
      print_assumption_line(x$assumptions)
      tf_blank()
    }
    tf_section("Result")
    tf_line(paste0("Raw n = ", format_sample_value(x$n)))
    if (!is.null(x$n_adjusted)) {
      tf_line(paste0("Adjusted n = ", format_sample_size_value(x$n_adjusted)))
    }
    if (!is.null(x$n_per_group)) {
      tf_line(paste0("Adjusted per-group n = ", format_sample_size_value(x$n_per_group)))
    }
    if (!is.null(x$n_total)) {
      tf_line(paste0("Adjusted total n = ", format_sample_size_value(x$n_total)))
    }
    if (!is.null(x$required_events)) {
      tf_line(paste0("Required events = ", format_sample_count(x$required_events)))
    }
    tf_blank()
    tf_section("Report")
    tf_line(report(x))
    print_diagnostics_section(x$diagnostics)
    invisible(x)
  })
}

# Compact "Diagnostics" section shared by print.sample_size() and
# print.summary.sample_size(); a no-op when diagnostics are absent (e.g. for
# endpoints/designs other than one-sample binary precision).
print_diagnostics_section <- function(d) {
  if (is.null(d)) return(invisible(NULL))
  tf_blank()
  tf_section("Diagnostics")
  tf_field("Anticipated prevalence", format_sample_value(d$anticipated_prevalence))
  tf_field("Confidence level", format_sample_value(d$confidence_level))
  tf_field("CI method", d$ci_method)
  tf_field("Precision criterion", d$precision_criterion)
  tf_field("Complete-case n", format_sample_count(d$complete_case_n))
  tf_field("Adjusted (recruitment) n", format_sample_count(d$adjusted_n))
  tf_field("Expected events (adjusted n)", format_sample_value(d$expected_events))
  tf_field("P(zero events)", format_sample_value(d$probability_zero_events))
  tf_field("Achieved maximum half-width", format_sample_value(d$achieved_maximum_half_width))
  if (!is.null(d$approximation_warning)) {
    tf_line(paste0("Note: ", d$approximation_warning))
  }
  invisible(NULL)
}

#' Summarize a sample size object
#'
#' @param object A sample-size object.
#' @param ... Unused.
#' @return A compact summary list.
#' @export
summary.sample_size <- function(object, ...) {
  out <- list(
    endpoint = object$endpoint,
    design = object$design,
    objective = object$objective,
    method = object$method,
    n = object$n,
    n_adjusted = object$n_adjusted,
    n_per_group = object$n_per_group,
    n_total = object$n_total,
    required_events = object$required_events,
    assumptions = object$assumptions,
    report = report(object),
    diagnostics = object$diagnostics
  )
  class(out) <- "summary.sample_size"
  out
}

#' Print a summary for a sample size object
#'
#' @param x A `summary.sample_size` object.
#' @param ... Unused.
#' @return The object invisibly.
#' @export
print.summary.sample_size <- function(x, ...) {
  tf_with_cli_colors({
    tf_title("sample size summary")
    print_summary_field("Endpoint", x$endpoint)
    print_summary_field("Design", x$design)
    print_summary_field("Objective", x$objective)
    print_summary_field("Method", x$method)
    print_summary_field("Raw n", x$n, formatter = format_sample_value)
    print_summary_field("Adjusted n", x$n_adjusted, formatter = format_sample_size_value)
    print_summary_field("Adjusted per-group n", x$n_per_group, formatter = format_sample_size_value)
    print_summary_field("Adjusted total n", x$n_total, formatter = format_sample_size_value)
    print_summary_field("Required events", x$required_events, formatter = format_sample_count)
    if (!is.null(x$report)) {
      tf_blank()
      tf_section("Report")
      tf_line(x$report)
    }
    print_diagnostics_section(x$diagnostics)
    invisible(x)
  })
}

#' @noRd
#' @export
as_tibble.sample_size <- function(x, ...) {
  out <- .bt(
    endpoint = x$endpoint,
    design = x$design,
    objective = x$objective,
    method = x$method,
    n = x$n,
    n_adjusted = format_sample_size_value(x$n_adjusted),
    n_per_group = format_sample_size_value(x$n_per_group),
    n_total = format_sample_size_value(x$n_total),
    required_events = format_sample_size_value(x$required_events),
    report = x$report
  )
  d <- x$diagnostics
  if (!is.null(d)) {
    out <- .as_bt(cbind(out, .bt(
      anticipated_prevalence = d$anticipated_prevalence,
      confidence_level = d$confidence_level,
      ci_method = d$ci_method,
      precision_criterion = d$precision_criterion,
      complete_case_n = d$complete_case_n,
      adjusted_n = d$adjusted_n,
      expected_events = d$expected_events,
      expected_non_events = d$expected_non_events,
      probability_zero_events = d$probability_zero_events,
      probability_at_least_one_event = d$probability_at_least_one_event,
      achieved_maximum_half_width = d$achieved_maximum_half_width,
      dropout = d$dropout
    )))
  }
  out
}

#' Return a sample size report
#'
#' @param x A sample-size object.
#' @param ... Unused.
#' @return A length-one character vector.
#' @export
report.sample_size <- function(x, ...) {
  x$report
}

sample_size_result <- function(endpoint, design, objective, method, n, n_adjusted, n_per_group = NULL, n_total = NULL, required_events = NULL, assumptions = NULL, report, plot_data = NULL, curve_data = NULL, diagnostics = NULL) {
  x <- list(
    endpoint = endpoint,
    design = design,
    objective = objective,
    method = method,
    n = n,
    n_adjusted = n_adjusted,
    n_per_group = n_per_group,
    n_total = n_total,
    required_events = required_events,
    assumptions = assumptions,
    report = report,
    plot_data = plot_data,
    curve_data = curve_data,
    diagnostics = diagnostics
  )
  # Planning inputs of the calling samplesize*() function, kept so plot() can draw a power curve for any design.
  pf <- parent.frame()
  x$alpha <- get0("alpha", pf, mode = "numeric", inherits = FALSE)
  x$power <- get0("power", pf, mode = "numeric", inherits = FALSE)
  # The calling function and its inputs let plot(vary = ) rerun the design under other assumptions;
  # `powerfun` (exact power against the total n) is set where the design has one (noncentral t).
  x$fun <- sys.function(-1); x$inputs <- ssinputs(x$fun, pf)
  x$powerfun <- attr(n, "powerfun"); attr(x$n, "powerfun") <- NULL
  class(x) <- c("sample_size", "testflow_sample_size")
  x
}

# Arguments given to (or defaulted in) a samplesize*() call, as their current values; missing ones are left out.
ssinputs <- function(fn, env) {
  nm <- setdiff(names(formals(fn)), "...")
  nm <- nm[!vapply(nm, function(v) eval(call("missing", as.name(v)), env) && is.symbol(formals(fn)[[v]]), TRUE)]
  mget(nm, envir = env)
}

sample_size_plot_frame <- function(labels, raw_values, adjusted_values) {
  labels <- as.character(labels)
  if (length(labels) != length(raw_values) || length(labels) != length(adjusted_values)) {
    stop("Plot labels and values must have the same length.", call. = FALSE)
  }
  raw_df <- .bt(label = labels, type = "raw", value = as.numeric(raw_values), display = vapply(raw_values, format_sample_size_value, character(1)))
  adj_df <- .bt(label = labels, type = "adjusted", value = as.numeric(adjusted_values), display = vapply(adjusted_values, format_sample_size_value, character(1)))
  .rbind_rows(raw_df, adj_df)
}

sample_size_curve_bundle <- function(n_target, power_target, n_grid, power_grid, unit_label, power_label) {
  .bt(
    n = as.numeric(n_grid),
    power = as.numeric(power_grid),
    target_n = as.numeric(n_target),
    power_target = as.numeric(power_target),
    n_label = unit_label,
    power_label = power_label
  )
}

sample_size_curve_grid <- function(n_target, n_min = NULL, n_max = NULL, length.out = 20) {
  n_target <- max(1, as.numeric(n_target))
  if (is.null(n_min)) n_min <- max(2, floor(n_target * 0.25))
  if (is.null(n_max)) n_max <- max(10, ceiling(n_target * 1.35))
  seq.int(n_min, n_max, length.out = length.out)
}

sample_size_curve_power_continuous_superiority <- function(n, delta, sd, allocation = 1, alpha = 0.05) {
  sample_size_validate_positive(delta, "delta")
  sample_size_validate_positive(sd, "sd")
  sample_size_validate_positive(allocation, "allocation")
  z_alpha <- stats::qnorm(1 - alpha / 2)
  stats::pnorm(sqrt(n * allocation / (allocation + 1)) * delta / sd - z_alpha)
}

sample_size_curve_power_continuous_paired_superiority <- function(n, delta, sd, alpha = 0.05) {
  sample_size_validate_positive(delta, "delta")
  sample_size_validate_positive(sd, "sd")
  z_alpha <- stats::qnorm(1 - alpha / 2)
  stats::pnorm(sqrt(n) * delta / sd - z_alpha)
}

sample_size_repeated_effective_sd <- function(sd, n_time, correlation) {
  sample_size_validate_positive(sd, "sd")
  sample_size_validate_positive(n_time, "n_time")
  sample_size_validate_probability(correlation, "correlation", allow_zero = TRUE)
  sd * sqrt((1 + (n_time - 1) * correlation) / n_time)
}

sample_size_assumptions <- function(...) {
  items <- list(...)
  items <- Filter(function(.x) !is.null(.x) && length(.x) > 0 && !all(is.na(.x)) && nzchar(as.character(.x)[1]), items)
  .imap_dfr(items, function(x, idx) {
    assumption_check(
      name = paste0("Planning note ", idx),
      status = "assumed",
      message = as.character(x)
    )
  })
}

sample_size_validate_probability <- function(x, name, allow_zero = FALSE) {
  if (length(x) != 1 || is.na(x) || !is.numeric(x)) {
    stop("`", name, "` must be a single numeric value.", call. = FALSE)
  }
  if (allow_zero) {
    if (x < 0 || x >= 1) stop("`", name, "` must be in [0, 1).", call. = FALSE)
  } else if (x <= 0 || x >= 1) {
    stop("`", name, "` must be in (0, 1).", call. = FALSE)
  }
  invisible(x)
}

sample_size_validate_positive <- function(x, name) {
  if (length(x) != 1 || is.na(x) || !is.numeric(x) || x <= 0) {
    stop("`", name, "` must be a single positive numeric value.", call. = FALSE)
  }
  invisible(x)
}

sample_size_validate_dropout <- function(dropout) {
  sample_size_validate_probability(dropout, "dropout", allow_zero = TRUE)
}

sample_size_validate_nonnegative <- function(x, name) {
  if (length(x) != 1 || is.na(x) || !is.numeric(x) || x < 0) {
    stop("`", name, "` must be a single non-negative numeric value.", call. = FALSE)
  }
  invisible(x)
}

sample_size_validate_max_n <- function(max_n) {
  if (length(max_n) != 1 || is.na(max_n) || !is.numeric(max_n) || !is.finite(max_n) || max_n < 2 || !isTRUE(all.equal(max_n, round(max_n)))) {
    stop("`max_n` must be a single finite integer greater than or equal to 2.", call. = FALSE)
  }
  as.integer(round(max_n))
}

# Smallest n (per arm A for parallel designs, pairs for paired designs) whose
# noncentral-t power reaches the target; `n_start` is the normal-approximation n.
sample_size_continuous_t <- function(n_start, objective, d, sd, alpha, power, paired = FALSE, allocation = 1, margin = NULL) {
  power_at <- function(n) {
    if (paired) {
      se <- sd / sqrt(n)
      df <- n - 1
    } else {
      n_b <- allocation * n
      se <- sd * sqrt(1 / n + 1 / n_b)
      df <- n + n_b - 2
    }
    if (df < 1) return(0)
    switch(
      objective,
      superiority = {
        tc <- stats::qt(1 - alpha / 2, df)
        stats::pt(tc, df, ncp = d / se, lower.tail = FALSE) + stats::pt(-tc, df, ncp = d / se)
      },
      noninferiority = stats::pt(stats::qt(1 - alpha, df), df, ncp = d / se, lower.tail = FALSE),
      equivalence = {
        tc <- stats::qt(1 - alpha, df)
        max(0, stats::pt(-tc, df, ncp = (d - margin) / se) - stats::pt(tc, df, ncp = (d + margin) / se))
      }
    )
  }
  n <- max(2, floor(n_start))
  while (power_at(n) < power) n <- n + 1
  # Exact power against the total size (pairs, or both arms), kept for the power curve.
  structure(n, powerfun = function(N) power_at(if (paired) N else N / (1 + allocation)))
}

sample_size_round <- function(x) {
  as.integer(ceiling(x))
}

sample_size_apply_dropout <- function(n, dropout) {
  sample_size_round(n / (1 - dropout))
}

sample_size_events_to_n <- function(events_total, survival_a = NULL, survival_b = NULL, accrual_duration = NULL, follow_up = NULL) {
  if (is.null(survival_a) || is.null(survival_b)) {
    return(events_total)
  }
  sample_size_validate_probability(survival_a, "survival_a", allow_zero = TRUE)
  sample_size_validate_probability(survival_b, "survival_b", allow_zero = TRUE)
  if (!is.null(accrual_duration) || !is.null(follow_up)) {
    if (is.null(accrual_duration) || is.null(follow_up)) {
      stop("Provide both `accrual_duration` and `follow_up` for accrual-adjusted event probabilities.", call. = FALSE)
    }
    sample_size_validate_positive(accrual_duration, "accrual_duration")
    sample_size_validate_positive(follow_up, "follow_up")
    total_duration <- accrual_duration + follow_up
    q_a <- sample_size_uniform_accrual_event_prob(survival_a, accrual_duration, total_duration)
    q_b <- sample_size_uniform_accrual_event_prob(survival_b, accrual_duration, total_duration)
  } else {
    q_a <- 1 - survival_a
    q_b <- 1 - survival_b
  }
  2 * events_total / (q_a + q_b)
}

# Uniform-accrual event probability (spec Sec. 24.2): survival is assumed
# exponential with hazard implied by the planned survival probability at the
# total study duration (accrual + follow-up), and enrollment is uniform over
# the accrual window. Reduces to `1 - survival` as accrual_duration -> 0
# (instantaneous accrual), matching the non-accrual conversion above.
sample_size_uniform_accrual_event_prob <- function(survival, accrual_duration, total_duration) {
  if (survival <= 0) return(1)
  lambda <- -log(survival) / total_duration
  1 - (exp(-lambda * (total_duration - accrual_duration)) - exp(-lambda * total_duration)) / (lambda * accrual_duration)
}

format_sample_count <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (length(x) > 1) {
    rendered <- vapply(unname(x), function(z) format_sample_count(z), character(1))
    if (is.null(names(x)) || all(!nzchar(names(x)))) {
      return(paste(rendered, collapse = ", "))
    }
    vals <- paste0(names(x), "=", rendered)
    return(paste(vals, collapse = ", "))
  }
  if (is.na(x)) return(NA_character_)
  format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
}

format_sample_value <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (length(x) > 1) {
    rendered <- vapply(unname(x), function(z) format_sample_value(z), character(1))
    if (is.null(names(x)) || all(!nzchar(names(x)))) {
      return(paste(rendered, collapse = ", "))
    }
    vals <- paste0(names(x), "=", rendered)
    return(paste(vals, collapse = ", "))
  }
  if (is.na(x)) return(NA_character_)
  formatC(x, format = "f", digits = 3)
}

format_sample_size_value <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (length(x) > 1) {
    rendered <- vapply(unname(x), function(z) format_sample_size_value(z), character(1))
    if (is.null(names(x)) || all(!nzchar(names(x)))) {
      return(paste(rendered, collapse = ", "))
    }
    vals <- paste0(names(x), "=", rendered)
    return(paste(vals, collapse = ", "))
  }
  if (is.na(x)) return(NA_character_)
  format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
}
