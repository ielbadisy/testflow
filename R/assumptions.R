#' Check normality
#' @param data A data frame.
#' @param vars Numeric columns as character names.
#' @param group Optional grouping column name.
#' @param alpha Significance level.
#' @noRd
check_normality <- function(data, vars, group = NULL, alpha = 0.05) {
  run_one <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 3) {
      return(assumption_check("Normality", "not enough data", "Too few observations for a useful normality check.", method = "Shapiro-Wilk", statistic = NA_real_, p_value = NA_real_))
    }
    if (n > 5000) {
      warning("Shapiro-Wilk is not recommended for n > 5000.", call. = FALSE)
      return(assumption_check("Normality", "not checked", "Shapiro-Wilk is not recommended for n > 5000.", method = "Shapiro-Wilk", statistic = NA_real_, p_value = NA_real_))
    }
    test <- stats::shapiro.test(x)
    assumption_check(
      "Normality",
      ifelse(test$p.value >= alpha, "acceptable", "not acceptable"),
      ifelse(test$p.value >= alpha, "Approximate normality looks reasonable.", "Normality may be violated."),
      method = "Shapiro-Wilk",
      statistic = unname(test$statistic),
      p_value = test$p.value
    )
  }

  vars <- as.character(vars)
  if (is.null(group)) {
    purrr::map_dfr(vars, function(v) {
      dplyr::mutate(run_one(data[[v]]), name = paste0("Normality: ", v), variable = v)
    })
  } else {
    purrr::map_dfr(vars, function(v) {
      split(data[[v]], data[[group]]) |>
        purrr::imap_dfr(function(x, g) {
          dplyr::mutate(run_one(x), name = paste0("Normality: ", v, " (", g, ")"), variable = v, group = g)
        })
    })
  }
}

#' Check symmetry about the median (Cabilio-Masaro test)
#' @param x Numeric vector of deviations from the reference value.
#' @param alpha Significance level.
#' @noRd
check_symmetry <- function(x, alpha = 0.05) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 3 || stats::sd(x) == 0) {
    return(assumption_check("Symmetry of deviations", "not enough data", "Too few observations, or no variation, for a useful symmetry check.", method = "Cabilio-Masaro test"))
  }
  stat <- sqrt(n) * (mean(x) - stats::median(x)) / (stats::sd(x) * sqrt(pi / 2 - 1))
  p_value <- 2 * (1 - stats::pnorm(abs(stat)))
  assumption_check(
    "Symmetry of deviations",
    ifelse(p_value >= alpha, "acceptable", "not acceptable"),
    ifelse(p_value >= alpha, "The deviations from the reference value look approximately symmetric.", "The deviations from the reference value show significant asymmetry."),
    method = "Cabilio-Masaro test",
    statistic = stat,
    p_value = p_value
  )
}

#' Check monotonicity
#' @noRd
check_monotonicity <- function(x, y, alpha = 0.05) {
  x <- x[stats::complete.cases(x, y)]
  y <- y[stats::complete.cases(x, y)]
  if (length(x) < 3) {
    return(assumption_check("Monotonic relationship", "not checked", "Not enough complete observations to assess monotonicity."))
  }
  sp <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  assumption_check(
    "Monotonic relationship",
    ifelse(sp$p.value >= alpha, "acceptable", "warning"),
    ifelse(sp$p.value >= alpha, "A monotonic association appears plausible.", "Relationship may be non-monotonic."),
    method = "Spearman correlation",
    statistic = unname(sp$statistic),
    p_value = sp$p.value
  )
}

#' Check linearity via a RESET-style quadratic-term F-test
#' @param x,y Numeric vectors.
#' @param alpha Significance level.
#' @noRd
check_linearity <- function(x, y, alpha = 0.05) {
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 4) {
    return(assumption_check("Linearity", "not checked", "Not enough complete observations to assess linearity.", method = "Quadratic-term F-test"))
  }
  fit_linear <- stats::lm(y ~ x)
  fit_quad <- stats::lm(y ~ x + I(x^2))
  cmp <- stats::anova(fit_linear, fit_quad)
  stat <- cmp$F[2]
  p_value <- cmp[["Pr(>F)"]][2]
  assumption_check(
    "Linearity",
    ifelse(p_value >= alpha, "acceptable", "warning"),
    ifelse(p_value >= alpha, "Adding a quadratic term does not significantly improve fit; a linear relation looks reasonable.", "Adding a quadratic term significantly improves fit, suggesting curvature; inspect a scatterplot."),
    method = "Quadratic-term F-test (Ramsey RESET, power = 2)",
    statistic = stat,
    p_value = p_value
  )
}

#' Check outliers
#' @noRd
check_outliers <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 4) {
    return(assumption_check("Extreme outliers", "not checked", "Too few observations for a useful outlier screen."))
  }
  q1 <- stats::quantile(x, 0.25, names = FALSE)
  q3 <- stats::quantile(x, 0.75, names = FALSE)
  i <- q3 - q1
  n <- sum(x < q1 - 1.5 * i | x > q3 + 1.5 * i)
  assumption_check("Extreme outliers", ifelse(n == 0, "acceptable", "warning"), paste0(n, " potential outlier(s) flagged by IQR."), details = paste0("IQR rule, n = ", n))
}

#' Check expected counts
#'
#' Uses the same "any expected cell below `threshold`" rule as the
#' chi-square/Fisher recommendation in `test_categorical()` (see
#' `fisher_threshold`), so this panel explains the recommendation actually
#' made rather than applying a different (Cochran 80%-of-cells) convention
#' that could disagree with it.
#' @noRd
check_expected_counts <- function(tab, threshold = 5) {
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  expected <- as.vector(chi$expected)
  ok <- all(expected >= threshold)
  assumption_check(
    "Expected cell counts",
    ifelse(ok, "acceptable", "not acceptable"),
    ifelse(ok, "Chi-square approximation is reasonable.", "Some expected counts are below the threshold; Fisher's exact test is used instead."),
    method = "Pearson chi-square approximation",
    details = paste0("Min expected = ", format(min(expected), digits = 3), "; threshold = ", format(threshold, digits = 3))
  )
}

#' Check sphericity via Mauchly's test
#' @param wide_matrix A numeric matrix, one row per subject, one column per
#'   repeated-measures condition (no missing values).
#' @param alpha Significance level.
#' @noRd
check_sphericity <- function(wide_matrix, alpha = 0.05) {
  k <- ncol(wide_matrix)
  if (k < 3) {
    return(assumption_check(
      "Sphericity", "not applicable",
      "Sphericity is not defined with only two repeated-measures conditions; a single within-subject contrast has no covariance structure to test.",
      method = "Mauchly's test"
    ))
  }
  mlm_fit <- stats::lm(wide_matrix ~ 1)
  mt <- stats::mauchly.test(mlm_fit, X = ~1)
  assumption_check(
    "Sphericity",
    ifelse(mt$p.value >= alpha, "acceptable", "warning"),
    ifelse(mt$p.value >= alpha, "Sphericity looks reasonable.", "Sphericity may be violated; interpret the repeated-measures ANOVA with caution, or prefer the Friedman test."),
    method = "Mauchly's test",
    statistic = unname(mt$statistic),
    p_value = unname(mt$p.value)
  )
}

#' Check homogeneity of variance
#' @noRd
check_variance_homogeneity <- function(data, outcome, group, alpha = 0.05) {
  formula <- stats::as.formula(paste(outcome, "~", group))
  lt <- car::leveneTest(formula, data = data)
  p <- lt[["Pr(>F)"]][1]
  stat <- lt[["F value"]][1]
  tibble::tibble(
    name = "Variance homogeneity",
    method = "Levene test",
    message = ifelse(p >= alpha, "Variance homogeneity looks reasonable.", "Variance homogeneity may be violated."),
    statistic = stat,
    df1 = lt[["Df"]][1],
    df2 = lt[["Df"]][2],
    p = p,
    status = ifelse(p >= alpha, "acceptable", "not acceptable")
  )
}

#' Check variance for two groups
#' @noRd
check_variance_two_groups <- function(data, outcome, group, alpha = 0.05) {
  formula <- stats::as.formula(paste(outcome, "~", group))
  test <- stats::var.test(formula, data = data)
  tibble::tibble(
    name = "Variance ratio check",
    message = ifelse(test$p.value >= alpha, "Variance ratio looks reasonable.", "Variance ratio looks concerning."),
    statistic = unname(test$statistic),
    df1 = unname(test$parameter[1]),
    df2 = unname(test$parameter[2]),
    p = test$p.value,
    conf.low = test$conf.int[1],
    conf.high = test$conf.int[2],
    status = ifelse(test$p.value >= alpha, "acceptable", "not acceptable")
  )
}

#' Check Bartlett homogeneity
#' @noRd
check_bartlett <- function(data, outcome, group, alpha = 0.05) {
  test <- stats::bartlett.test(stats::as.formula(paste(outcome, "~", group)), data = data)
  tibble::tibble(
    name = "Bartlett test",
    method = "Bartlett test",
    message = ifelse(test$p.value >= alpha, "Variance homogeneity looks reasonable.", "Variance homogeneity may be violated."),
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p = test$p.value,
    status = ifelse(test$p.value >= alpha, "acceptable", "not acceptable")
  )
}
#' Build a standardized assumption check
#' @noRd
assumption_check <- function(name, status, message, method = NA_character_, statistic = NA_real_, p_value = NA_real_, details = NA_character_, ...) {
  tibble::tibble(
    name = name,
    status = status,
    message = message,
    method = method,
    statistic = statistic,
    p_value = p_value,
    details = details,
    ...
  )
}

#' Combine assumption checks
#' @noRd
assumption_checks <- function(...) {
  dplyr::bind_rows(...)
}

#' Format assumption checks
#' @noRd
format_assumptions <- function(assumptions) {
  if (is.null(assumptions)) {
    return(tibble::tibble())
  }
  if (is.data.frame(assumptions)) {
    return(tibble::as_tibble(assumptions))
  }
  if (is.list(assumptions)) {
    out <- purrr::imap_dfr(assumptions, function(x, nm) {
      if (is.data.frame(x)) {
        x$name <- x$name %||% nm
        return(x)
      }
      assumption_check(name = nm, status = "not checked", message = as.character(x))
    })
    return(out)
  }
  tibble::tibble()
}

#' Check independence note
#' @noRd
check_independence_note <- function(message = "Assumed from study design.") {
  assumption_check("Independence of observations", "assumed", message)
}

#' Check normality
