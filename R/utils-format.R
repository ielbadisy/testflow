#' Format a p-value
#' @keywords internal
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
#' @keywords internal
format_stat <- function(x, digits = 2) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

#' Classify Cohen's d magnitude
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
magnitude_cramers_v <- function(v) {
  dplyr::case_when(
    is.na(v) ~ NA_character_,
    v < 0.1 ~ "negligible",
    v < 0.3 ~ "small",
    v < 0.5 ~ "moderate",
    TRUE ~ "large"
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
