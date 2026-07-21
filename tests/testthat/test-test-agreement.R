test_that("test_agreement matches the Fleiss-Cohen-Everitt (1969) kappa SE formula exactly", {
  set.seed(1)
  n <- 100
  rater1 <- sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  rater2 <- rater1
  flip <- sample(seq_len(n), 20)
  rater2[flip] <- sample(c("A", "B", "C"), 20, replace = TRUE)
  dat <- tibble::tibble(rater1 = rater1, rater2 = rater2)

  x <- test_agreement(dat, rater1, rater2)
  expect_s3_class(x, "testflow_agreement")

  tab <- table(rater1, rater2)
  n_tot <- sum(tab)
  p <- tab / n_tot
  pi_row <- rowSums(p); pi_col <- colSums(p)
  po <- sum(diag(p)); pe <- sum(pi_row * pi_col)
  kappa_ref <- (po - pe) / (1 - pe)
  expect_equal(x$effect_size$estimate[1], kappa_ref)

  k <- nrow(p)
  term1 <- sum(sapply(seq_len(k), function(i) p[i, i] * (1 - (pi_row[i] + pi_col[i]) * (1 - kappa_ref))^2))
  term2 <- (1 - kappa_ref)^2 * sum(sapply(seq_len(k), function(i) sum(sapply(seq_len(k), function(j) if (i == j) 0 else p[i, j] * (pi_col[i] + pi_row[j])^2))))
  term3 <- (kappa_ref - pe * (1 - kappa_ref))^2
  se_ref <- sqrt((term1 + term2 - term3) / (n_tot * (1 - pe)^2))
  ci_ref <- kappa_ref + c(-1, 1) * stats::qnorm(0.975) * se_ref
  expect_equal(c(x$primary_test$conf.low, x$primary_test$conf.high), ci_ref)

  # The z-test of H0: kappa = 0 must use the NULL-hypothesis SE (Fleiss,
  # 1981), evaluated at kappa = 0 - a different, generally smaller-magnitude
  # SE than se_ref above (which is appropriate for the CI, not this test).
  # Reusing se_ref here was a real bug: it understates the null SE whenever
  # kappa is materially above 0, inflating z and anti-conservatively
  # shrinking the p-value.
  se0_ref <- sqrt(pe + pe^2 - sum(pi_row * pi_col * (pi_row + pi_col))) / ((1 - pe) * sqrt(n_tot))
  z_ref <- kappa_ref / se0_ref
  p_ref <- 2 * stats::pnorm(-abs(z_ref))
  expect_equal(x$primary_test$statistic, z_ref)
  expect_equal(x$primary_test$p.value, p_ref)
  expect_equal(x$effect_size$magnitude[1], "substantial")

  # se0_ref must differ from se_ref (otherwise this test can't distinguish
  # the fix from the original bug), and must be the larger of the two here
  # (kappa is well above chance, the typical case where the bug mattered).
  expect_true(se0_ref > se_ref)

  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
})

test_that("test_agreement's kappa hypothesis test matches irr::kappa2()'s null-hypothesis z-test", {
  skip_if_not_installed("irr")
  set.seed(2)
  n <- 200
  cats <- c("A", "B", "C")
  rater1 <- sample(cats, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  rater2 <- ifelse(stats::runif(n) < 0.75, rater1, sample(cats, n, replace = TRUE))
  dat <- data.frame(r1 = factor(rater1, levels = cats), r2 = factor(rater2, levels = cats))

  x <- test_agreement(dat, rater1 = r1, rater2 = r2)
  irr_res <- irr::kappa2(dat[, c("r1", "r2")])

  expect_equal(unname(x$effect_size$estimate[1]), unname(irr_res$value), tolerance = 1e-8)
  expect_equal(unname(x$primary_test$statistic), unname(irr_res$statistic), tolerance = 1e-6)
  expect_equal(unname(x$primary_test$p.value), unname(irr_res$p.value), tolerance = 1e-6)
})

test_that("test_agreement handles category sets that appear in only one rater's column", {
  dat <- tibble::tibble(
    rater1 = c("yes", "yes", "no", "no", "no"),
    rater2 = c("yes", "no", "no", "no", "maybe")
  )
  x <- test_agreement(dat, rater1, rater2)
  expect_s3_class(x, "testflow_agreement")
  expect_equal(nrow(x$alternative_tests$agreement_table), 3)
})

test_that("test_icc matches independently computed Shrout-Fleiss/McGraw-Wong ANOVA formulas exactly", {
  set.seed(1)
  n <- 20; k <- 3
  true_score <- rnorm(n, 50, 10)
  ratings <- sapply(seq_len(k), function(j) true_score + rnorm(n, 0, 5))
  colnames(ratings) <- paste0("r", seq_len(k))
  dat <- as.data.frame(ratings)

  x <- test_icc(dat, c(r1, r2, r3))
  expect_s3_class(x, "testflow_icc")

  long <- data.frame(subject = factor(rep(seq_len(n), k)), rater = factor(rep(seq_len(k), each = n)), value = as.vector(ratings))
  fit1 <- aov(value ~ subject, data = long)
  tab1 <- anova(fit1)
  bms <- tab1["subject", "Mean Sq"]; wms <- tab1["Residuals", "Mean Sq"]
  icc1_ref <- (bms - wms) / (bms + (k - 1) * wms)

  fit2 <- aov(value ~ subject + rater, data = long)
  tab2 <- anova(fit2)
  bms2 <- tab2["subject", "Mean Sq"]; jms <- tab2["rater", "Mean Sq"]; ems <- tab2["Residuals", "Mean Sq"]
  icc3_ref <- (bms2 - ems) / (bms2 + (k - 1) * ems)
  icc2_ref <- (bms2 - ems) / (bms2 + (k - 1) * ems + (k / n) * (jms - ems))

  tbl <- x$alternative_tests$icc_table
  expect_equal(tbl$estimate[1], icc1_ref)
  expect_equal(tbl$estimate[2], icc2_ref)
  expect_equal(tbl$estimate[3], icc3_ref)
  expect_equal(x$effect_size$estimate[1], icc2_ref)

  # F-test for ICC(2,1) (H0: ICC = 0) collapses to BMS/EMS on (n-1),(n-1)(k-1) df -
  # the same F-statistic as ICC(3,1)'s test, since a=0,b=1 at the null.
  f_ref <- bms2 / ems
  expect_equal(tbl$statistic[2], f_ref)
  expect_equal(tbl$p.value[2], stats::pf(f_ref, n - 1, (n - 1) * (k - 1), lower.tail = FALSE))
  expect_equal(x$primary_test$statistic, f_ref)

  expect_equal(x$effect_size$magnitude[1], "good")
  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
})

test_that("test_icc's ICC(2,1) confidence interval uses the estimated ICC, not the null value", {
  set.seed(1)
  n <- 20; k <- 3
  true_score <- rnorm(n, 50, 10)
  ratings <- sapply(seq_len(k), function(j) true_score + rnorm(n, 0, 5))
  colnames(ratings) <- paste0("r", seq_len(k))
  dat <- as.data.frame(ratings)
  x <- test_icc(dat, c(r1, r2, r3))
  tbl <- x$alternative_tests$icc_table
  # The CI must be centered near the point estimate, not near 0 (which would
  # indicate the null-hypothesis coefficients leaked into the CI computation).
  expect_true(tbl$conf.low[2] < tbl$estimate[2])
  expect_true(tbl$conf.high[2] > tbl$estimate[2])
  expect_true(tbl$conf.low[2] > 0.5)
})
