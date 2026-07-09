#' Format a p-value
#' @noRd
format_p <- function(p, digits = 3) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(
      p < 10^-digits,
      paste0("<", formatC(10^-digits, format = "f", digits = digits)),
      formatC(p, format = "f", digits = digits)
    )
  )
}

#' Format a statistic
#' @noRd
format_stat <- function(x, digits = 2) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

#' Classify Cohen's d magnitude
#' @noRd
magnitude_cohens_d <- function(d) {
  ad <- abs(d)
  dplyr::case_when(
    is.na(ad) ~ NA_character_,
    ad < 0.2 ~ "negligible",
    ad < 0.5 ~ "small",
    ad < 0.8 ~ "moderate",
    TRUE ~ "large"
  )
}

#' Classify eta squared magnitude
#' @noRd
magnitude_eta2 <- function(eta2) {
  dplyr::case_when(
    is.na(eta2) ~ NA_character_,
    eta2 < 0.01 ~ "negligible",
    eta2 < 0.06 ~ "small",
    eta2 < 0.14 ~ "moderate",
    TRUE ~ "large"
  )
}

#' Classify Cramer's V magnitude
#' @noRd
magnitude_cramers_v <- function(v) {
  dplyr::case_when(
    is.na(v) ~ NA_character_,
    v < 0.1 ~ "negligible",
    v < 0.3 ~ "small",
    v < 0.5 ~ "moderate",
    TRUE ~ "large"
  )
}

#' Classify AUC magnitude (Hosmer, Lemeshow & Sturdivant, 2013, p. 177)
#' @noRd
magnitude_auc <- function(auc) {
  a <- ifelse(is.na(auc), NA_real_, pmax(auc, 1 - auc))
  dplyr::case_when(
    is.na(a) ~ NA_character_,
    a < 0.7 ~ "no better than chance / poor",
    a < 0.8 ~ "acceptable",
    a < 0.9 ~ "excellent",
    TRUE ~ "outstanding"
  )
}

#' Classify Cohen's kappa magnitude (Landis & Koch, 1977)
#' @noRd
magnitude_kappa <- function(kappa) {
  dplyr::case_when(
    is.na(kappa) ~ NA_character_,
    kappa < 0 ~ "poor",
    kappa < 0.21 ~ "slight",
    kappa < 0.41 ~ "fair",
    kappa < 0.61 ~ "moderate",
    kappa < 0.81 ~ "substantial",
    TRUE ~ "almost perfect"
  )
}

#' Classify ICC magnitude (Koo & Li, 2016)
#' @noRd
magnitude_icc <- function(icc) {
  dplyr::case_when(
    is.na(icc) ~ NA_character_,
    icc < 0.5 ~ "poor",
    icc < 0.75 ~ "moderate",
    icc < 0.9 ~ "good",
    TRUE ~ "excellent"
  )
}

significance_decision <- function(p, alpha = 0.05) {
  if (is.na(p)) {
    "No decision because the p-value is unavailable."
  } else if (p < alpha) {
    "Statistically significant at the chosen alpha level."
  } else {
    "Not statistically significant at the chosen alpha level."
  }
}
