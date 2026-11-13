test_that("testdiagnostic matches hand-computed sensitivity/specificity/LR and CIs exactly", {
  tp <- 45; fp <- 10; fn <- 8; tn <- 90
  dat <- .bt(
    test = c(rep("positive", tp + fp), rep("negative", fn + tn)),
    reference = c(rep("positive", tp), rep("negative", fp), rep("positive", fn), rep("negative", tn))
  )
  x <- testdiagnostic(dat, test, reference)
  expect_s3_class(x, "testflow_diagnostic")

  tbl <- x$alternative_tests$diagnostic_table
  expect_equal(tbl$estimate[tbl$metric == "Sensitivity"], tp / (tp + fn))
  expect_equal(tbl$estimate[tbl$metric == "Specificity"], tn / (tn + fp))
  expect_equal(tbl$estimate[tbl$metric == "Positive predictive value"], tp / (tp + fp))
  expect_equal(tbl$estimate[tbl$metric == "Negative predictive value"], tn / (tn + fn))
  expect_equal(tbl$estimate[tbl$metric == "Accuracy"], (tp + tn) / (tp + fp + fn + tn))

  sens_ci_ref <- as.numeric(stats::binom.test(tp, tp + fn)$conf.int)
  expect_equal(c(tbl$conf.low[tbl$metric == "Sensitivity"], tbl$conf.high[tbl$metric == "Sensitivity"]), sens_ci_ref)

  lrp_ref <- (tp / (tp + fn)) / (1 - tn / (tn + fp))
  se_ref <- sqrt(1 / tp - 1 / (tp + fn) + 1 / fp - 1 / (fp + tn))
  lrp_ci_ref <- exp(log(lrp_ref) + c(-1, 1) * stats::qnorm(0.975) * se_ref)
  expect_equal(tbl$estimate[tbl$metric == "Positive likelihood ratio"], lrp_ref)
  expect_equal(c(tbl$conf.low[tbl$metric == "Positive likelihood ratio"], tbl$conf.high[tbl$metric == "Positive likelihood ratio"]), lrp_ci_ref)

  n <- tp + fp + fn + tn
  nir <- max((tp + fn) / n, (fp + tn) / n)
  acc_ref <- stats::binom.test(tp + tn, n, p = nir, alternative = "greater")
  expect_equal(x$primary_test$statistic, unname(acc_ref$statistic))
  expect_equal(x$primary_test$p.value, acc_ref$p.value)

  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
})

test_that("testdiagnostic requires exactly two levels for test and reference", {
  dat <- .bt(test = sample(c("a", "b", "c"), 30, replace = TRUE), reference = sample(c("pos", "neg"), 30, replace = TRUE))
  expect_error(testdiagnostic(dat, test, reference), "exactly two")
})

test_that("testroc matches wilcox.test-derived AUC, Hanley-McNeil CI, and a trapezoidal cross-check", {
  set.seed(1)
  n1 <- 50; n2 <- 60
  dat <- .bt(
    marker = c(rnorm(n2, 2, 1), rnorm(n1, 0, 1)),
    disease = c(rep("yes", n2), rep("no", n1))
  )
  x <- testroc(dat, marker, disease)
  expect_s3_class(x, "testflow_roc")

  wt <- stats::wilcox.test(dat$marker[dat$disease == "yes"], dat$marker[dat$disease == "no"])
  auc_ref <- unname(wt$statistic) / (n1 * n2)
  expect_equal(x$effect_size$estimate[1], auc_ref)

  q1 <- auc_ref / (2 - auc_ref)
  q2 <- 2 * auc_ref^2 / (1 + auc_ref)
  se_ref <- sqrt((auc_ref * (1 - auc_ref) + (n2 - 1) * (q1 - auc_ref^2) + (n1 - 1) * (q2 - auc_ref^2)) / (n1 * n2))
  ci_ref <- auc_ref + c(-1, 1) * stats::qnorm(0.975) * se_ref
  expect_equal(c(x$primary_test$conf.low, x$primary_test$conf.high), ci_ref)

  z_ref <- (auc_ref - 0.5) / se_ref
  p_ref <- 2 * stats::pnorm(-abs(z_ref))
  expect_equal(x$primary_test$p.value, p_ref)

  rc <- x$alternative_tests$roc_curve
  fpr <- 1 - rc$specificity
  auc_trap <- sum(diff(fpr) * (head(rc$sensitivity, -1) + tail(rc$sensitivity, -1)) / 2)
  expect_equal(auc_trap, auc_ref, tolerance = 1e-10)

  expect_equal(x$effect_size$magnitude[1], "outstanding")
  expect_s3_class(plot(x), "ggplot")
  expect_type(report(x), "character")
})

test_that("testroc's Youden's-J threshold is the ROC point maximizing sensitivity + specificity - 1", {
  set.seed(2)
  n1 <- 40; n2 <- 40
  dat <- .bt(
    marker = c(rnorm(n2, 1, 1), rnorm(n1, 0, 1)),
    disease = c(rep("yes", n2), rep("no", n1))
  )
  x <- testroc(dat, marker, disease)
  rc <- x$alternative_tests$roc_curve
  j <- rc$sensitivity + rc$specificity - 1
  opt <- x$alternative_tests$optimal_threshold
  expect_equal(opt$youden_j, max(j))
})

test_that("testroc method = 'delong' gives the DeLong (1988) SE, default stays Hanley-McNeil", {
  set.seed(11)
  pos <- c(round(rnorm(40, 1), 1), 0.5, 0.5)
  neg <- c(round(rnorm(60, 0), 1), 0.5)
  dat <- data.frame(marker = c(pos, neg), disease = rep(c("yes", "no"), c(length(pos), length(neg))))

  # DeLong structural components by brute force (ties count 1/2)
  psi <- outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b))
  v10 <- rowMeans(psi); v01 <- colMeans(psi)
  se_delong <- sqrt(stats::var(v10) / length(pos) + stats::var(v01) / length(neg))
  z <- stats::qnorm(0.975)

  xd <- testroc(dat, marker, disease, plot = FALSE, method = "delong")
  auc <- mean(psi)
  expect_equal(xd$effect_size$estimate, auc, tolerance = 1e-12)
  expect_equal(xd$primary_test$conf.low, auc - z * se_delong, tolerance = 1e-10)
  expect_equal(xd$primary_test$conf.high, auc + z * se_delong, tolerance = 1e-10)
  expect_match(xd$primary_test$method, "DeLong")

  xh <- testroc(dat, marker, disease, plot = FALSE)
  xh2 <- testroc(dat, marker, disease, plot = FALSE, method = "hanley")
  expect_equal(xh$primary_test$conf.low, xh2$primary_test$conf.low)
  expect_false(isTRUE(all.equal(xh$primary_test$conf.low, xd$primary_test$conf.low)))
  expect_error(testroc(dat, marker, disease, plot = FALSE, method = "bogus"))
})
