#' Sample size planning for common trial designs
#'
#' @description
#' `sample_size()` dispatches to endpoint-specific planning helpers for
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
sample_size <- function(endpoint = c("continuous", "binary", "survival", "ordinal"), design = c("parallel", "paired", "repeated"), objective = c("superiority", "noninferiority", "equivalence"), ...) {
  endpoint <- match.arg(endpoint)
  design <- match.arg(design)
  objective <- match.arg(objective)
  switch(
    endpoint,
    continuous = sample_size_continuous(design = design, objective = objective, ...),
    binary = sample_size_binary(design = design, objective = objective, ...),
    survival = sample_size_survival(design = design, objective = objective, ...),
    ordinal = sample_size_ordinal(design = design, objective = objective, ...)
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
sample_size_continuous <- function(
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
  correlation = 0.5
) {
  design <- match.arg(design)
  objective <- match.arg(objective)
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
      n_main <- sample_size_round(n_raw)
      n_adj <- sample_size_apply_dropout(n_main, dropout)
      return(sample_size_result(
        endpoint = "continuous",
        design = design_label,
        objective = objective,
        method = if (repeated_design) "repeated-measures normal approximation" else "paired normal approximation",
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
      n_main <- sample_size_round(n_raw)
      n_adj <- sample_size_apply_dropout(n_main, dropout)
      return(sample_size_result(
        endpoint = "continuous",
        design = design_label,
        objective = objective,
        method = if (repeated_design) "repeated-measures non-inferiority normal approximation" else "paired non-inferiority normal approximation",
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
      n_main <- sample_size_round(n_raw)
      n_adj <- sample_size_apply_dropout(n_main, dropout)
      return(sample_size_result(
        endpoint = "continuous",
        design = design_label,
        objective = objective,
        method = if (repeated_design) "repeated-measures equivalence approximation" else "paired equivalence approximation",
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
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "continuous",
      design = "parallel",
      objective = objective,
      method = "parallel normal approximation",
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
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "continuous",
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
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "continuous",
      design = "parallel",
      objective = objective,
      method = "parallel equivalence approximation",
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
sample_size_binary <- function(
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
sample_size_survival <- function(
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
#'   target two-one-sided-test power exactly; the preferred method) or
#'   `normal_approx` (closed-form approximation).
#' @return A `sample_size` object.
#' @details
#' On the log scale, with \eqn{\theta_0=\log(GMR)}, \eqn{\theta_L=\log(L)},
#' \eqn{\theta_U=\log(U)}, and \eqn{d_{BE}=\min(\theta_0-\theta_L,\
#' \theta_U-\theta_0)} (must be positive: the anticipated GMR must lie
#' strictly inside the bioequivalence bounds):
#'
#' `method = "iterative_tost"` searches for the smallest \eqn{n} such that
#' the exact two-one-sided-test power is at least the target:
#' \deqn{Power(n) = \Phi\left(\frac{\theta_U-\theta_0}{SE(n)}-z_{1-\alpha}\right)
#' + \Phi\left(\frac{\theta_0-\theta_L}{SE(n)}-z_{1-\alpha}\right) - 1}
#' with \eqn{SE(n) = \sqrt{2\sigma_w^2/n}} for a crossover design (total
#' \eqn{n}) or \eqn{SE(n_A) = \sigma_b\sqrt{(r+1)/(rn_A)}} for a parallel
#' design (per-arm \eqn{n_A}, \eqn{n_B=rn_A}). Coefficients of variation are
#' converted to log-scale standard deviations via
#' \eqn{\sigma = \sqrt{\log(1+CV^2)}}.
#'
#' `method = "normal_approx"` uses the closed-form approximation
#' \deqn{n_{total} \approx \frac{2\sigma_w^2(z_{1-\beta}+z_{1-\alpha})^2}{d_{BE}^2}}
#' (crossover) or, per arm,
#' \deqn{n_A \approx \frac{(r+1)\sigma_b^2(z_{1-\beta}+z_{1-\alpha})^2}{rd_{BE}^2}}
#' (parallel), switching to \eqn{z_{1-\beta/2}} in place of \eqn{z_{1-\beta}}
#' when \eqn{GMR=1} exactly, matching the analogous special case in
#' [sample_size_continuous()]. `iterative_tost` is preferred because this
#' approximation can meaningfully differ from the exact TOST power away
#' from \eqn{GMR=1}.
#' @references
#' Julious SA. *Sample Sizes for Clinical Trials*. Chapman & Hall/CRC; 2010.
#'
#' Phillips KF. Power of the two one-sided tests procedure in
#' bioequivalence. *Journal of Pharmacokinetics and Biopharmaceutics*.
#' 1990;18(2):137-144.
#' @export
sample_size_bioequivalence <- function(
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
  } else {
    sample_size_validate_positive(cv_between, "cv_between")
    sigma <- sqrt(log(1 + cv_between^2))
    se_fn <- function(n) sigma * sqrt((allocation + 1) / (allocation * n))
  }

  tost_power <- function(n) {
    se <- se_fn(n)
    stats::pnorm((theta_u - theta0) / se - z_alpha_one) + stats::pnorm((theta0 - theta_l) / se - z_alpha_one) - 1
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
    n_main <- sample_size_round(n_raw)
    n_adj <- sample_size_apply_dropout(n_main, dropout)
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
sample_size_ordinal <- function(
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

#' Adjust a sample size for dropout
#'
#' @param n Evaluability sample size.
#' @param dropout Dropout proportion.
#' @return The dropout-adjusted sample size.
#' @export
sample_size_adjust_dropout <- function(n, dropout = 0) {
  sample_size_validate_dropout(dropout)
  sample_size_round(n / (1 - dropout))
}

#' Adjust an individually randomized sample size for cluster randomization
#'
#' @param n_ind Individually randomized sample size (e.g. from
#'   [sample_size_continuous()] or [sample_size_binary()]).
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
sample_size_cluster_adjust <- function(n_ind, m, rho, cv_m = NULL) {
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

#' Sample size for confidence-interval precision
#'
#' @param endpoint Endpoint family: `continuous` or `binary`.
#' @param design `one_sample`, `two_sample`, or (binary only) `odds_ratio`.
#' @param width Desired confidence-interval half-width, on the scale of the
#'   estimate (or of the log odds ratio for `design = "odds_ratio"`).
#' @param sd Standard deviation (continuous endpoint).
#' @param p Anticipated proportion (binary, `one_sample`).
#' @param p1 Anticipated proportion in arm A (binary, `two_sample` or
#'   `odds_ratio`).
#' @param p2 Anticipated proportion in arm B (binary, `two_sample` or
#'   `odds_ratio`).
#' @param alpha Type I error rate (two-sided CI coverage is `1 - alpha`).
#' @param allocation Allocation ratio `n_B / n_A`. Only used for
#'   `endpoint = "continuous"`, `design = "two_sample"` and
#'   `design = "odds_ratio"`; the binary `two_sample` formula is
#'   equal-allocation only.
#' @param dropout Expected dropout proportion.
#' @param conservative For `endpoint = "binary"`, `design = "one_sample"`
#'   only: use the worst-case `p = 0.5` instead of the anticipated `p`.
#' @return A `sample_size` object.
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
#' Binary, one proportion:
#' \deqn{n = \frac{z_{1-\alpha/2}^2p(1-p)}{w^2}}
#' (or \eqn{n = z_{1-\alpha/2}^2/(4w^2)} for the conservative \eqn{p=0.5} case)
#'
#' Binary, two independent proportions (equal allocation):
#' \deqn{n_{per\ group} = \frac{z_{1-\alpha/2}^2\left[p_A(1-p_A)+p_B(1-p_B)\right]}{w^2}}
#'
#' Binary, log-odds-ratio half-width:
#' \deqn{n_A = \frac{z_{1-\alpha/2}^2}{w^2}
#' \left[\frac1{p_A}+\frac1{1-p_A}+\frac1{rp_B}+\frac1{r(1-p_B)}\right]}
#' @references
#' Julious SA. *Sample Sizes for Clinical Trials*. Chapman & Hall/CRC; 2010.
#' @export
sample_size_precision <- function(
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
  conservative = FALSE
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
    if (conservative) {
      n_raw <- z_alpha^2 / (4 * width^2)
    } else {
      sample_size_validate_probability(p, "p")
      n_raw <- z_alpha^2 * p * (1 - p) / width^2
    }
    n_adj <- sample_size_apply_dropout(sample_size_round(n_raw), dropout)
    return(sample_size_result(
      endpoint = endpoint, design = design, objective = "precision",
      method = if (conservative) "one-proportion CI half-width (conservative p = 0.5)" else "one-proportion CI half-width",
      n = n_raw, n_adjusted = n_adj, n_total = n_adj,
      assumptions = sample_size_assumptions(paste0("Target CI half-width = ", format_sample_value(width), ".")),
      report = paste0("One-proportion precision sample size was calculated with width = ", format_sample_value(width), "; required n = ", format_sample_count(n_adj), "."),
      plot_data = sample_size_plot_frame("Required n", n_raw, n_adj)
    ))
  }

  sample_size_validate_probability(p1, "p1")
  sample_size_validate_probability(p2, "p2")
  if (identical(design, "two_sample")) {
    n_raw <- z_alpha^2 * (p1 * (1 - p1) + p2 * (1 - p2)) / width^2
    n_adj <- sample_size_apply_dropout(sample_size_round(n_raw), dropout)
    return(sample_size_result(
      endpoint = endpoint, design = design, objective = "precision",
      method = "two-proportion CI half-width (equal allocation)",
      n = n_raw, n_adjusted = c(A = n_adj, B = n_adj), n_per_group = c(A = n_adj, B = n_adj), n_total = 2 * n_adj,
      assumptions = sample_size_assumptions(paste0("Target CI half-width for the risk difference = ", format_sample_value(width), "."), "Equal allocation is assumed; this formula does not generalize to unequal allocation."),
      report = paste0("Two-proportion precision sample size was calculated with p1 = ", format_sample_value(p1), ", p2 = ", format_sample_value(p2), ", width = ", format_sample_value(width), "; required per-group size = ", format_sample_count(n_adj), "."),
      plot_data = sample_size_plot_frame(c("A", "B"), c(n_raw, n_raw), c(n_adj, n_adj))
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

#' Plot a sample size object
#'
#' @param x A sample-size object returned by `sample_size()`.
#' @param type Plot style: `both`, `curve`, or `summary`.
#' @param ... Unused.
#' @return A `ggplot2` object, or a combined patchwork object when `type =
#'   "both"`.
#' @export
plot.sample_size <- function(x, type = c("both", "curve", "summary"), ...) {
  type <- match.arg(type)
  summary_plot <- NULL
  curve_plot <- NULL

  if (type != "curve" && !is.null(x$plot_data) && nrow(x$plot_data)) {
    summary_plot <- ggplot2::ggplot(x$plot_data, ggplot2::aes(x = .data$label, y = .data$value, fill = .data$type)) +
      ggplot2::geom_col(width = 0.7, color = "white", position = ggplot2::position_dodge(width = 0.75)) +
      ggplot2::geom_text(ggplot2::aes(label = .data$display), vjust = -0.35, size = 3.6, position = ggplot2::position_dodge(width = 0.75)) +
      ggplot2::labs(
        title = paste0("Sample size planning: ", x$endpoint, " / ", x$design),
        subtitle = paste0(x$objective, " using ", x$method),
        x = NULL,
        y = "Sample size"
      ) +
      ggplot2::scale_fill_manual(values = c(raw = "#4C78A8", adjusted = "#F58518")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none")
  }

  if (type != "summary" && !is.null(x$curve_data) && is.data.frame(x$curve_data) && nrow(x$curve_data)) {
    label_target_n <- x$curve_data$target_n[1]
    label_y <- max(x$curve_data$power, na.rm = TRUE) * 0.18
    label_text <- paste0("optimal sample size\nn = ", format_sample_count(label_target_n))
    label_df <- tibble::tibble(x = label_target_n, y = label_y, label = label_text)
    curve_plot <- ggplot2::ggplot(x$curve_data, ggplot2::aes(x = .data$n, y = .data$power)) +
      ggplot2::geom_line(color = "#4C78A8", linewidth = 1.0) +
      ggplot2::geom_point(color = "#1f1f1f", size = 1.9) +
      ggplot2::geom_vline(xintercept = x$curve_data$target_n[1], linetype = "dotted", color = "#5B5BD6", linewidth = 0.9) +
      ggplot2::geom_label(data = label_df, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), inherit.aes = FALSE, hjust = 0, vjust = 0, linewidth = 0, fill = "white", alpha = 0.85, color = "#5B5BD6") +
      ggplot2::geom_hline(yintercept = x$curve_data$power_target[1], linetype = "dashed", color = "#F58518", linewidth = 0.7) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(formatC(100 * x, format = "f", digits = 0), "%"), limits = c(0, 1.05)) +
      ggplot2::labs(
        title = paste0("Sample size planning: ", x$endpoint, " / ", x$design),
        subtitle = paste0(x$objective, " using ", x$method, " | ", x$curve_data$power_label[1]),
        x = x$curve_data$n_label[1],
        y = "test power = 1 - beta"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  }

  if (type == "summary") return(summary_plot)
  if (type == "curve") return(curve_plot %||% summary_plot)
  if (!is.null(curve_plot) && !is.null(summary_plot)) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(patchwork::wrap_plots(curve_plot, summary_plot, ncol = 1))
    }
    return(curve_plot)
  }
  curve_plot %||% summary_plot
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
    invisible(x)
  })
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
    report = report(object)
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
    invisible(x)
  })
}

#' Convert a sample size object to a tibble
#'
#' @param x A sample-size object.
#' @param ... Unused.
#' @return A one-row tibble summarizing the planning result.
#' @export
as_tibble.sample_size <- function(x, ...) {
  tibble::tibble(
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

sample_size_result <- function(endpoint, design, objective, method, n, n_adjusted, n_per_group = NULL, n_total = NULL, required_events = NULL, assumptions = NULL, report, plot_data = NULL, curve_data = NULL) {
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
    curve_data = curve_data
  )
  class(x) <- c("sample_size", "testflow_sample_size")
  x
}

sample_size_plot_frame <- function(labels, raw_values, adjusted_values) {
  labels <- as.character(labels)
  if (length(labels) != length(raw_values) || length(labels) != length(adjusted_values)) {
    stop("Plot labels and values must have the same length.", call. = FALSE)
  }
  raw_df <- tibble::tibble(label = labels, type = "raw", value = as.numeric(raw_values), display = vapply(raw_values, format_sample_size_value, character(1)))
  adj_df <- tibble::tibble(label = labels, type = "adjusted", value = as.numeric(adjusted_values), display = vapply(adjusted_values, format_sample_size_value, character(1)))
  dplyr::bind_rows(raw_df, adj_df)
}

sample_size_curve_bundle <- function(n_target, power_target, n_grid, power_grid, unit_label, power_label) {
  tibble::tibble(
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
  items <- purrr::keep(items, ~ !is.null(.x) && length(.x) > 0 && !all(is.na(.x)) && nzchar(as.character(.x)[1]))
  purrr::imap_dfr(items, function(x, idx) {
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
