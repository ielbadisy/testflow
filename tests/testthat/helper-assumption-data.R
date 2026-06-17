make_assumption_data <- function(n = 80) {
  set.seed(1)
  sex <- rep(c("female", "male"), length.out = n)
  treatment <- rep(c("control", "treated"), length.out = n)
  id <- seq_len(n)
  sbp_baseline <- rnorm(n, 135, 12)
  sbp_3m <- sbp_baseline + ifelse(treatment == "treated", -5, 0) + rnorm(n, 0, 8)
  sbp_6m <- sbp_3m + rnorm(n, 0, 6)
  ldl <- rnorm(n, 120, 20)
  crp <- rlnorm(n, 1, 0.4)
  controlled_baseline <- sample(c("yes", "no"), n, replace = TRUE)
  controlled_3m <- sample(c("yes", "no"), n, replace = TRUE)
  controlled_6m <- sample(c("yes", "no"), n, replace = TRUE)
  age <- round(rnorm(n, 55, 10))
  tibble::tibble(
    id = id,
    age = age,
    sex = sex,
    treatment = treatment,
    sbp_baseline = sbp_baseline,
    sbp_3m = sbp_3m,
    sbp_6m = sbp_6m,
    ldl = ldl,
    crp = crp,
    controlled_baseline = controlled_baseline,
    controlled_3m = controlled_3m,
    controlled_6m = controlled_6m
  )
}
