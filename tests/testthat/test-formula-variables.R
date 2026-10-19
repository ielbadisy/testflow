test_that("a formula stored in a variable is accepted", {
  set.seed(1)
  df <- data.frame(
    y = rnorm(60), x = rnorm(60),
    g2 = rep(c("a", "b"), 30), g3 = rep(c("a", "b", "c"), 20),
    c2 = rep(c("u", "v"), each = 30), after = rnorm(60)
  )
  f2 <- y ~ g2
  expect_equal(test_two_groups(f2, df, plot = FALSE)$primary_test$p.value,
               test_two_groups(y ~ g2, df, plot = FALSE)$primary_test$p.value)
  expect_equal(test_two_groups(df, f2, plot = FALSE)$primary_test$p.value,
               test_two_groups(y ~ g2, df, plot = FALSE)$primary_test$p.value)
  f3 <- y ~ g3
  expect_equal(test_groups(f3, df, plot = FALSE)$primary_test$p.value,
               test_groups(y ~ g3, df, plot = FALSE)$primary_test$p.value)
  fc <- y ~ x
  expect_equal(test_correlation(fc, df, plot = FALSE)$primary_test$p.value,
               test_correlation(y ~ x, df, plot = FALSE)$primary_test$p.value)
  fk <- g2 ~ c2
  expect_equal(test_categorical(fk, df, plot = FALSE)$primary_test$p.value,
               test_categorical(g2 ~ c2, df, plot = FALSE)$primary_test$p.value)
  fp <- after ~ y
  expect_equal(test_paired(fp, df, plot = FALSE)$primary_test$p.value,
               test_paired(after ~ y, df, plot = FALSE)$primary_test$p.value)
  ff <- y ~ g2 * c2
  expect_equal(test_factorial(ff, df, plot = FALSE)$primary_test$p.value,
               test_factorial(y ~ g2 * c2, df, plot = FALSE)$primary_test$p.value)
  fs <- ~ y + g2 | c2
  expect_equal(sumtab(fs, df), sumtab(~ y + g2 | c2, df))
})

test_that("formulas built with as.formula() inline and inside functions work", {
  set.seed(2)
  df <- data.frame(y = rnorm(40), g = rep(c("a", "b"), 20))
  r1 <- test_two_groups(stats::as.formula("y ~ g"), df, plot = FALSE)
  wrap <- function(outcome, data) {
    f <- stats::as.formula(paste(outcome, "~ g"))
    test_two_groups(f, data, plot = FALSE)
  }
  expect_equal(wrap("y", df)$primary_test$p.value, r1$primary_test$p.value)
})

test_that("bare column names still work", {
  set.seed(3)
  df <- data.frame(y = rnorm(40), g = rep(c("a", "b"), 20))
  expect_equal(test_two_groups(df, y, g, plot = FALSE)$primary_test$p.value,
               test_two_groups(y ~ g, df, plot = FALSE)$primary_test$p.value)
})
