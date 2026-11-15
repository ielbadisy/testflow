test_that("testdiagnostic(ci = 'wilson') matches binom::binom.confint(method = 'wilson')", {
  skip_if_not_installed("binom")
  d <- data.frame(
    test = rep(c("neg", "pos", "neg", "pos"), c(80, 12, 9, 41)),
    ref = rep(c("absent", "absent", "present", "present"), c(80, 12, 9, 41))
  )
  w <- testdiagnostic(d, test, ref, ci = "wilson", plot = FALSE)
  tab <- w$alternative_tests$diagnostic_table
  counts <- list(c(41, 50), c(80, 92), c(41, 53), c(80, 89), c(121, 142))
  ref <- do.call(rbind, lapply(counts, function(k) binom::binom.confint(k[1], k[2], method = "wilson")[, c("lower", "upper")]))
  expect_equal(tab$conf.low[1:5], ref$lower)
  expect_equal(tab$conf.high[1:5], ref$upper)

  e <- testdiagnostic(d, test, ref, plot = FALSE)
  expect_equal(e$alternative_tests$diagnostic_table$conf.low[1], binom.test(41, 50)$conf.int[1])
  expect_equal(e$alternative_tests$diagnostic_table$estimate, tab$estimate)
})
