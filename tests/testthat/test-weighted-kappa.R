test_that("weighted kappa matches irr::kappa2 (estimate, z) and psych::cohen.kappa (CI)", {
  skip_if_not_installed("irr")
  skip_if_not_installed("psych")
  env <- new.env()
  utils::data("anxiety", package = "irr", envir = env)
  d <- data.frame(a = env$anxiety$rater1, b = env$anxiety$rater2)
  k <- length(unique(c(d$a, d$b)))
  for (w in c("linear", "quadratic")) {
    x <- testagreement(d, a, b, weights = w, plot = FALSE)
    ref <- irr::kappa2(d, weight = if (w == "linear") "equal" else "squared")
    expect_equal(x$effect_size$estimate[1], ref$value)
    expect_equal(x$primary_test$statistic, ref$statistic)
    wm <- if (w == "linear") 1 - abs(outer(1:k, 1:k, "-")) / (k - 1) else 1 - (outer(1:k, 1:k, "-") / (k - 1))^2
    ps <- suppressWarnings(psych::cohen.kappa(d, w = wm))
    expect_equal(c(x$primary_test$conf.low, x$primary_test$conf.high),
                 unname(ps$confid["weighted kappa", c("lower", "upper")]), tolerance = 1e-8)
    expect_match(x$effect_size$name[1], "Weighted kappa")
  }
})

test_that("weighted kappa follows factor level order and default stays unweighted", {
  lv <- c("none", "mild", "severe")
  r1 <- factor(c("none", "mild", "severe", "mild", "none", "severe", "mild", "none"), levels = lv)
  r2 <- factor(c("none", "severe", "severe", "mild", "mild", "mild", "mild", "none"), levels = lv)
  d <- data.frame(r1 = r1, r2 = r2)
  x <- testagreement(d, r1, r2, weights = "linear", plot = FALSE)
  ref <- irr::kappa2(data.frame(as.integer(r1), as.integer(r2)), weight = "equal")
  expect_equal(x$effect_size$estimate[1], ref$value)
  u <- testagreement(d, r1, r2, plot = FALSE)
  expect_equal(u$effect_size$name[1], "Cohen's kappa")
  expect_equal(u$effect_size$estimate[1], irr::kappa2(d)$value)
})
