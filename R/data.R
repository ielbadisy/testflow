#' Simulate a small cardiovascular teaching dataset
#' @param n Number of rows.
#' @param seed Random seed.
#' @return A tibble with example numeric and categorical variables.
#' @export
make_cardio_data <- function(n = 180, seed = 2026) {
  set.seed(seed)
  sex <- sample(c("female", "male"), n, replace = TRUE)
  treatment <- sample(c("usual care", "lifestyle", "medication"), n, replace = TRUE)
  smoker <- sample(c("no", "yes"), n, replace = TRUE, prob = c(0.72, 0.28))
  diabetes <- sample(c("no", "yes"), n, replace = TRUE, prob = c(0.78, 0.22))
  age <- round(stats::rnorm(n, 58, 11))
  sbp_baseline <- round(stats::rnorm(n, 146, 15) + ifelse(diabetes == "yes", 5, 0), 1)
  treatment_drop <- dplyr::case_when(treatment == "usual care" ~ 5, treatment == "lifestyle" ~ 10, TRUE ~ 14)
  sbp_3m <- round(sbp_baseline - treatment_drop + ifelse(sex == "male", 3, 0) + stats::rnorm(n, 0, 10), 1)
  sbp_6m <- round(sbp_3m - treatment_drop / 3 + stats::rnorm(n, 0, 9), 1)
  ldl <- round(stats::rnorm(n, 3.2, 0.8) + ifelse(smoker == "yes", 0.25, 0), 2)
  crp <- round(stats::rgamma(n, shape = 2, rate = 0.7), 2)
  adherence <- pmin(100, pmax(20, round(stats::rnorm(n, 78, 14))))
  response <- ifelse(sbp_3m < sbp_baseline - 10, "responder", "non-responder")
  tibble::tibble(
    id = seq_len(n),
    sex = sex,
    treatment = treatment,
    smoker = smoker,
    diabetes = diabetes,
    age = age,
    sbp_baseline = sbp_baseline,
    sbp_3m = sbp_3m,
    sbp_6m = sbp_6m,
    ldl = ldl,
    crp = crp,
    adherence = adherence,
    response = response,
    controlled_baseline = ifelse(sbp_baseline < 140, "yes", "no"),
    controlled_3m = ifelse(sbp_3m < 140, "yes", "no"),
    controlled_6m = ifelse(sbp_6m < 140, "yes", "no")
  )
}
