#' Check normality
#' @param data A data frame.
#' @param vars Numeric columns as character names.
#' @param group Optional grouping column name.
#' @param alpha Significance level.
#' @keywords internal
check_normality <- function(data, vars, group = NULL, alpha = 0.05) {
  run_one <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 3) {
      return(tibble::tibble(statistic = NA_real_, p = NA_real_, status = "not enough data"))
    }
    if (n > 5000) {
      warning("Shapiro-Wilk is not recommended for n > 5000.", call. = FALSE)
      return(tibble::tibble(statistic = NA_real_, p = NA_real_, status = "not checked"))
    }
    test <- stats::shapiro.test(x)
    tibble::tibble(
      statistic = unname(test$statistic),
      p = test$p.value,
      status = ifelse(test$p.value >= alpha, "acceptable", "not acceptable")
    )
  }

  vars <- as.character(vars)
  if (is.null(group)) {
    purrr::map_dfr(vars, function(v) {
      dplyr::bind_cols(
        tibble::tibble(
          name = "Normality",
          variable = v,
          group = NA_character_
        ),
        run_one(data[[v]])
      )
    })
  } else {
    purrr::map_dfr(vars, function(v) {
      split(data[[v]], data[[group]]) |>
        purrr::imap_dfr(function(x, g) {
          dplyr::bind_cols(
            tibble::tibble(
              name = "Normality",
              variable = v,
              group = as.character(g)
            ),
            run_one(x)
          )
        })
    })
  }
}

#' Check monotonicity
#' @keywords internal
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

#' Check outliers
#' @keywords internal
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
#' @keywords internal
check_expected_counts <- function(tab, threshold = 5) {
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
  expected <- as.vector(chi$expected)
  ok <- sum(expected >= threshold) / length(expected) >= 0.8 && all(expected >= 1)
  assumption_check(
    "Expected cell counts",
    ifelse(ok, "acceptable", "not acceptable"),
    ifelse(ok, "Chi-square approximation is reasonable.", "Some expected counts are too small for the chi-square approximation."),
    method = "Pearson chi-square approximation",
    details = paste0("Min expected = ", format(min(expected), digits = 3))
  )
}

#' Check sphericity or note
#' @keywords internal
check_sphericity_or_note <- function(...) {
  assumption_check("Sphericity", "not checked", "Sphericity is not checked here; use this as a teaching note unless a formal test is added.")
}

#' Check homogeneity of variance
#' @keywords internal
check_variance_homogeneity <- function(data, outcome, group, alpha = 0.05) {
  formula <- stats::as.formula(paste(outcome, "~", group))
  if (requireNamespace("car", quietly = TRUE)) {
    lt <- car::leveneTest(formula, data = data)
    p <- lt[["Pr(>F)"]][1]
    stat <- lt[["F value"]][1]
    return(tibble::tibble(
      method = "Levene test",
      statistic = stat,
      df1 = lt[["Df"]][1],
      df2 = lt[["Df"]][2],
      p = p,
      status = ifelse(p >= alpha, "acceptable", "not acceptable")
    ))
  }

  med <- stats::ave(data[[outcome]], data[[group]], FUN = stats::median, na.rm = TRUE)
  z <- abs(data[[outcome]] - med)
  fit <- stats::lm(z ~ data[[group]])
  aov_tab <- stats::anova(fit)
  p <- aov_tab[["Pr(>F)"]][1]
  tibble::tibble(
    method = "Median-centered Levene approximation",
    statistic = aov_tab[["F value"]][1],
    df1 = aov_tab[["Df"]][1],
    df2 = aov_tab[["Df"]][2],
    p = p,
    status = ifelse(p >= alpha, "acceptable", "not acceptable")
  )
}

#' Check variance for two groups
#' @keywords internal
check_variance_two_groups <- function(data, outcome, group, alpha = 0.05) {
  formula <- stats::as.formula(paste(outcome, "~", group))
  test <- stats::var.test(formula, data = data)
  tibble::tibble(
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
#' @keywords internal
check_bartlett <- function(data, outcome, group, alpha = 0.05) {
  test <- stats::bartlett.test(stats::as.formula(paste(outcome, "~", group)), data = data)
  tibble::tibble(
    method = "Bartlett test",
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p = test$p.value,
    status = ifelse(test$p.value >= alpha, "acceptable", "not acceptable")
  )
}
#' Build a standardized assumption check
#' @keywords internal
assumption_check <- function(name, status, message, method = NA_character_, statistic = NA_real_, p_value = NA_real_, details = NA_character_) {
  tibble::tibble(
    name = name,
    status = status,
    message = message,
    method = method,
    statistic = statistic,
    p_value = p_value,
    details = details
  )
}

#' Combine assumption checks
#' @keywords internal
assumption_checks <- function(...) {
  dplyr::bind_rows(...)
}

#' Format assumption checks
#' @keywords internal
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
#' @keywords internal
check_independence_note <- function(message = "Assumed from study design.") {
  assumption_check("Independence of observations", "assumed", message)
}

#' Check normality
