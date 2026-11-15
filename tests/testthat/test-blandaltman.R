test_that("testblandaltman matches the Bland-Altman formulas computed by hand", {
  set.seed(11)
  x <- rnorm(35, 100, 10)
  d <- data.frame(a = x + rnorm(35, 1.5, 4), b = x + rnorm(35, 0, 4))
  r <- testblandaltman(d, a, b, plot = FALSE)
  expect_s3_class(r, "testflow_blandaltman")

  dif <- d$a - d$b
  n <- length(dif); m <- mean(dif); s <- sd(dif)
  z <- qnorm(0.975); tq <- qt(0.975, n - 1)
  loa <- m + c(-1, 1) * z * s
  se_loa <- s * sqrt(1 / n + z^2 / (2 * (n - 1)))
  expect_equal(r$limits$estimate, c(m, loa))
  expect_equal(r$limits$conf.low, c(m - tq * s / sqrt(n), loa - tq * se_loa))
  expect_equal(r$limits$conf.high, c(m + tq * s / sqrt(n), loa + tq * se_loa))

  tt <- t.test(d$a, d$b, paired = TRUE)
  expect_equal(r$primary_test$p.value, tt$p.value)
  expect_equal(c(r$primary_test$conf.low, r$primary_test$conf.high), as.numeric(tt$conf.int))

  pb <- summary(lm(dif ~ I((d$a + d$b) / 2)))$coefficients[2, ]
  expect_equal(r$alternative_tests$proportional_bias$estimate, unname(pb[1]))
  expect_equal(r$alternative_tests$proportional_bias$p.value, unname(pb[4]))
})

test_that("testblandaltman alpha sets the LoA multiplier and the plot is a ggplot", {
  set.seed(3)
  d <- data.frame(a = rnorm(20), b = rnorm(20))
  r <- testblandaltman(d, a, b, alpha = 0.1)
  dif <- d$a - d$b
  expect_equal(r$limits$estimate[3], mean(dif) + qnorm(0.95) * sd(dif))
  expect_s3_class(r$plot, "ggplot")
  expect_match(report(r), "90% limits of agreement")
  expect_output(print(r), "Bias and limits of agreement")
})
