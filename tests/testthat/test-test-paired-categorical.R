# Tests for testpairedcategorical() (McNemar).

test_that("testpairedcategorical exposes the continuity correction and counts discordant pairs", {
  dat <- data.frame(
    before = factor(c(rep("yes", 30), rep("yes", 5), rep("no", 15), rep("no", 50)), levels = c("no", "yes")),
    after = factor(c(rep("yes", 30), rep("no", 5), rep("yes", 15), rep("no", 50)), levels = c("no", "yes"))
  )
  tab <- table(dat$before, dat$after)
  x <- testpairedcategorical(dat, before, after, plot = FALSE)
  expect_equal(x$primary_test$statistic, unname(stats::mcnemar.test(tab)$statistic))
  y <- testpairedcategorical(dat, before, after, correct = FALSE, plot = FALSE)
  expect_equal(y$primary_test$statistic, unname(stats::mcnemar.test(tab, correct = FALSE)$statistic))
  expect_equal(y$primary_test$statistic, (15 - 5)^2 / 20)
  disc <- x$assumptions[x$assumptions$name == "Discordant pairs", ]
  expect_match(disc$details, "discordant = 20", fixed = TRUE)
  expect_equal(x$effect_size$estimate, 20)
})
