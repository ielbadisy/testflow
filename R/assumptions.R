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
      dplyr::bind_cols(tibble::tibble(variable = v, group = NA_character_), run_one(data[[v]]))
    })
  } else {
    purrr::map_dfr(vars, function(v) {
      split(data[[v]], data[[group]]) |>
        purrr::imap_dfr(function(x, g) {
          dplyr::bind_cols(tibble::tibble(variable = v, group = as.character(g)), run_one(x))
        })
    })
  }
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
