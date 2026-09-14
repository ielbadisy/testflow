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
  as.character(cut(
    ad, breaks = c(-Inf, 0.2, 0.5, 0.8, Inf),
    labels = c("negligible", "small", "moderate", "large"), right = FALSE
  ))
}

#' Classify eta squared magnitude
#' @noRd
magnitude_eta2 <- function(eta2) {
  as.character(cut(
    eta2, breaks = c(-Inf, 0.01, 0.06, 0.14, Inf),
    labels = c("negligible", "small", "moderate", "large"), right = FALSE
  ))
}

#' Classify Cramer's V magnitude
#' @noRd
magnitude_cramers_v <- function(v) {
  as.character(cut(
    v, breaks = c(-Inf, 0.1, 0.3, 0.5, Inf),
    labels = c("negligible", "small", "moderate", "large"), right = FALSE
  ))
}

#' Classify AUC magnitude (Hosmer, Lemeshow & Sturdivant, 2013, p. 177)
#' @noRd
magnitude_auc <- function(auc) {
  a <- ifelse(is.na(auc), NA_real_, pmax(auc, 1 - auc))
  as.character(cut(
    a, breaks = c(-Inf, 0.7, 0.8, 0.9, Inf),
    labels = c("no better than chance / poor", "acceptable", "excellent", "outstanding"), right = FALSE
  ))
}

#' Classify Cohen's kappa magnitude (Landis & Koch, 1977)
#' @noRd
magnitude_kappa <- function(kappa) {
  as.character(cut(
    kappa, breaks = c(-Inf, 0, 0.21, 0.41, 0.61, 0.81, Inf),
    labels = c("poor", "slight", "fair", "moderate", "substantial", "almost perfect"), right = FALSE
  ))
}

#' Classify ICC magnitude (Koo & Li, 2016)
#' @noRd
magnitude_icc <- function(icc) {
  as.character(cut(
    icc, breaks = c(-Inf, 0.5, 0.75, 0.9, Inf),
    labels = c("poor", "moderate", "good", "excellent"), right = FALSE
  ))
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
