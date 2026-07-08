#' Sample size planning for common trial designs
#'
#' @description
#' `sample_size()` dispatches to endpoint-specific planning helpers. The
#' implementation covers the formulas selected in the Julious reference for
#' continuous, binary, survival, and ordinal planning problems.
#'
#' @param endpoint Endpoint family: `continuous`, `binary`, `survival`, or
#'   `ordinal`.
#' @param design Study design: `parallel`, `paired`, or `repeated`.
#' @param objective Planning objective: `superiority`, `noninferiority`,
#'   `equivalence`, or `precision` depending on endpoint.
#' @param ... Endpoint-specific arguments passed to the selected helper.
#' @return A `sample_size` object.
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
#' @return A `sample_size` object.
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
  n_time = 2
) {
  design <- match.arg(design)
  objective <- match.arg(objective)
  sample_size_validate_probability(alpha, "alpha")
  sample_size_validate_probability(power, "power")
  sample_size_validate_dropout(dropout)
  sample_size_validate_positive(allocation, "allocation")
  if (design == "repeated" && n_time > 2) {
    stop("Repeated-measures sample size is currently supported only when `n_time = 2`.", call. = FALSE)
  }

  z_power <- stats::qnorm(power)
  z_alpha_two <- stats::qnorm(1 - alpha / 2)
  z_alpha_one <- stats::qnorm(1 - alpha)
  paired_design <- design %in% c("paired", "repeated")

  if (paired_design) {
    scale_sd <- sd_diff %||% sd
    sample_size_validate_positive(scale_sd, "sd_diff")

    if (objective == "superiority") {
      sample_size_validate_positive(delta, "delta")
      n_raw <- ((z_alpha_two + z_power)^2 * scale_sd^2) / delta^2
      n_main <- sample_size_round(n_raw)
      n_adj <- sample_size_apply_dropout(n_main, dropout)
      return(sample_size_result(
        endpoint = "continuous",
        design = if (design == "repeated") "paired (two time points)" else "paired",
        objective = objective,
        method = "paired normal approximation",
        n = n_raw,
        n_adjusted = n_adj,
        n_total = n_adj,
        assumptions = sample_size_assumptions(
          "Paired observations are assumed.",
          paste0("Paired-difference SD = ", format_sample_value(scale_sd), ".")
        ),
        formula = "n = ((z_{1-alpha/2} + z_{1-beta})^2 * sd_diff^2) / delta^2",
        reference = SAMPLE_SIZE_REFERENCE,
        report = paste0(
          "Paired continuous superiority sample size was calculated with delta = ",
          format_sample_value(delta), ", sd_diff = ", format_sample_value(scale_sd),
          ", alpha = ", format_sample_value(alpha), ", power = ", format_sample_value(power),
          "; required pairs = ", format_sample_count(n_adj), "."
        ),
        plot_data = sample_size_plot_frame("Required pairs", n_raw, n_adj)
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
        design = if (design == "repeated") "paired (two time points)" else "paired",
        objective = objective,
        method = "paired non-inferiority normal approximation",
        n = n_raw,
        n_adjusted = n_adj,
        n_total = n_adj,
        assumptions = sample_size_assumptions(
          "Paired observations are assumed.",
          paste0("Distance to the non-inferiority boundary = ", format_sample_value(d_ni), ".")
        ),
        formula = "n = ((z_{1-alpha} + z_{1-beta})^2 * sd_diff^2) / (expected_difference + margin)^2",
        reference = SAMPLE_SIZE_REFERENCE,
        report = paste0(
          "Paired continuous non-inferiority sample size was calculated with expected_difference = ",
          format_sample_value(expected_difference), ", margin = ", format_sample_value(delta),
          ", sd_diff = ", format_sample_value(scale_sd), "; required pairs = ", format_sample_count(n_adj), "."
        ),
        plot_data = sample_size_plot_frame("Required pairs", n_raw, n_adj)
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
        formula_text <- "n = 2 * (z_{1-beta/2} + z_{1-alpha})^2 * sd_diff^2 / margin^2"
      } else {
        n_raw <- ((z_power + z_alpha_one)^2 * scale_sd^2) / d_eq^2
        formula_text <- "n = ((z_{1-beta} + z_{1-alpha})^2 * sd_diff^2) / (margin - abs(expected_difference))^2"
      }
      n_main <- sample_size_round(n_raw)
      n_adj <- sample_size_apply_dropout(n_main, dropout)
      return(sample_size_result(
        endpoint = "continuous",
        design = if (design == "repeated") "paired (two time points)" else "paired",
        objective = objective,
        method = "paired equivalence approximation",
        n = n_raw,
        n_adjusted = n_adj,
        n_total = n_adj,
        assumptions = sample_size_assumptions(
          "Paired observations are assumed.",
          paste0("Distance to the nearest equivalence boundary = ", format_sample_value(d_eq), ".")
        ),
        formula = formula_text,
        reference = SAMPLE_SIZE_REFERENCE,
        report = paste0(
          "Paired continuous equivalence sample size was calculated with expected_difference = ",
          format_sample_value(expected_difference), ", margin = ", format_sample_value(delta),
          ", sd_diff = ", format_sample_value(scale_sd), "; required pairs = ", format_sample_count(n_adj), "."
        ),
        plot_data = sample_size_plot_frame("Required pairs", n_raw, n_adj)
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
      formula = "n_A = ((r + 1) * sd^2 * (z_{1-beta} + z_{1-alpha/2})^2) / (r * delta^2)",
      reference = SAMPLE_SIZE_REFERENCE,
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
      formula = "n_A = ((r + 1) * sd^2 * (z_{1-beta} + z_{1-alpha})^2) / (r * (expected_difference + margin)^2)",
      reference = SAMPLE_SIZE_REFERENCE,
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
      formula_text <- "n_A = ((r + 1) * sd^2 * (z_{1-beta/2} + z_{1-alpha})^2) / (r * margin^2)"
    } else {
      n_a_raw <- ((allocation + 1) * sd^2 * (z_power + z_alpha_one)^2) / (allocation * d_eq^2)
      formula_text <- "n_A = ((r + 1) * sd^2 * (z_{1-beta} + z_{1-alpha})^2) / (r * (margin - abs(expected_difference))^2)"
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
      formula = formula_text,
      reference = SAMPLE_SIZE_REFERENCE,
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
#' @export
sample_size_binary <- function(
  design = c("parallel", "paired", "repeated"),
  objective = c("superiority", "noninferiority", "equivalence"),
  p1,
  p2,
  margin = NULL,
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
      formula = "n_total = ((z_{1-alpha/2} + z_{1-beta})^2 * (OR_D + 1)^2) / ((OR_D - 1)^2 * discordance_rate)",
      reference = SAMPLE_SIZE_REFERENCE,
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
    n_a_raw <- if (is.null(margin)) n_pooled else n_alt
    n_b_raw <- allocation * n_a_raw
    n_a <- sample_size_round(n_a_raw)
    n_b <- sample_size_round(n_b_raw)
    n_a_adj <- sample_size_apply_dropout(n_a, dropout)
    n_b_adj <- sample_size_apply_dropout(n_b, dropout)
    return(sample_size_result(
      endpoint = "binary",
      design = "parallel",
      objective = objective,
      method = if (is.null(margin)) "parallel pooled normal approximation" else "parallel anticipated-response normal approximation",
      n = n_a_raw,
      n_adjusted = c(A = n_a_adj, B = n_b_adj),
      n_per_group = c(A = n_a_adj, B = n_b_adj),
      assumptions = sample_size_assumptions(
        "Two independent groups are assumed.",
        paste0("Risk difference = ", format_sample_value(delta), ".")
      ),
      formula = if (is.null(margin)) {
        "n_A = [z_{1-alpha/2} * sqrt((1 + 1/r) * pbar * (1 - pbar)) + z_{1-beta} * sqrt(p1(1-p1) + p2(1-p2)/r)]^2 / (p1 - p2)^2"
      } else {
        "n_A = (z_{1-alpha/2} + z_{1-beta})^2 * [p1(1-p1) + p2(1-p2)/r] / (p1 - p2)^2"
      },
      reference = SAMPLE_SIZE_REFERENCE,
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
      formula = "n_A = (z_{1-alpha} + z_{1-beta})^2 * [p1(1-p1) + p2(1-p2)/r] / (p1 - p2 + margin)^2",
      reference = SAMPLE_SIZE_REFERENCE,
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
      n_a_raw <- ((z_power_eq + z_alpha_one)^2 * p1 * (1 - p1)) / (margin^2) * 2
      formula_text <- "n_per_group = 2 * (z_{1-beta/2} + z_{1-alpha})^2 * p(1-p) / margin^2"
    } else {
      n_a_raw <- ((z_power + z_alpha_one)^2 * (p1 * (1 - p1) + p2 * (1 - p2) / allocation)) / d_eq^2
      formula_text <- "n_A = (z_{1-alpha} + z_{1-beta})^2 * [p1(1-p1) + p2(1-p2)/r] / (margin - abs(p1 - p2))^2"
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
      formula = formula_text,
      reference = SAMPLE_SIZE_REFERENCE,
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
#' @return A `sample_size` object.
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
  method = c("exponential", "ph_only")
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
    n_total_raw <- sample_size_events_to_n(total_events, survival_a, survival_b)
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
        "Event probabilities are derived from the planned survival probabilities."
      ),
      formula = if (method == "exponential") {
        "D_total = 4 * (z_{1-alpha/2} + z_{1-beta})^2 / log(hr)^2"
      } else {
        "D_total = 2 * ((hr + 1)^2 * (z_{1-alpha/2} + z_{1-beta})^2) / (hr - 1)^2"
      },
      reference = SAMPLE_SIZE_REFERENCE,
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
    n_total_raw <- sample_size_events_to_n(total_events, survival_a, survival_b)
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
        paste0("Distance to the non-inferiority boundary = ", format_sample_value(d_ni), ".")
      ),
      formula = "D_total = 4 * (z_{1-alpha} + z_{1-beta})^2 / (log(margin_hr) - log(hr))^2",
      reference = SAMPLE_SIZE_REFERENCE,
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
    n_total_raw <- sample_size_events_to_n(total_events, survival_a, survival_b)
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
        paste0("Distance to the nearest equivalence boundary = ", format_sample_value(d_eq), ".")
      ),
      formula = "D_total = 4 * (z_{1-alpha} + z_{1-beta})^2 / d_EQ^2",
      reference = SAMPLE_SIZE_REFERENCE,
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
#' @export
sample_size_ordinal <- function(
  design = c("parallel"),
  objective = c("superiority"),
  p_superiority,
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
  sample_size_validate_probability(p_superiority, "p_superiority")
  if (p_superiority <= 0.5) {
    stop("p_superiority must exceed 0.5 for superiority planning.", call. = FALSE)
  }
  z_power <- stats::qnorm(power)
  z_alpha_two <- stats::qnorm(1 - alpha / 2)
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
    formula = "n_per_group = (z_{1-alpha/2} + z_{1-beta})^2 / [6 * (P - 0.5)^2]",
    reference = SAMPLE_SIZE_REFERENCE,
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

#' Plot a sample size object
#'
#' @param x A sample-size object returned by `sample_size()`.
#' @param ... Unused.
#' @return A `ggplot2` object.
#' @export
plot.sample_size <- function(x, ...) {
  if (is.null(x$plot_data) || !nrow(x$plot_data)) {
    return(NULL)
  }
  ggplot2::ggplot(x$plot_data, ggplot2::aes(x = .data$label, y = .data$value, fill = .data$type)) +
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
    formula = object$formula,
    reference = object$reference,
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
    formula = x$formula,
    reference = x$reference,
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

SAMPLE_SIZE_REFERENCE <- "Julious SA. Sample Sizes for Clinical Trials. Chapman & Hall/CRC; 2010."

sample_size_result <- function(endpoint, design, objective, method, n, n_adjusted, n_per_group = NULL, n_total = NULL, required_events = NULL, assumptions = NULL, formula, reference, report, plot_data = NULL) {
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
    formula = formula,
    reference = reference,
    report = report,
    plot_data = plot_data
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

sample_size_events_to_n <- function(events_total, survival_a = NULL, survival_b = NULL) {
  if (is.null(survival_a) || is.null(survival_b)) {
    return(events_total)
  }
  sample_size_validate_probability(survival_a, "survival_a", allow_zero = TRUE)
  sample_size_validate_probability(survival_b, "survival_b", allow_zero = TRUE)
  q_a <- 1 - survival_a
  q_b <- 1 - survival_b
  events_total / (q_a + q_b)
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
