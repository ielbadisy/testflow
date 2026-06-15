test_that("sumtab builds grouped numeric and categorical summaries", {
  dat <- tibble::tibble(
    age = c(10, 12, 20, 22),
    sex = factor(c("female", "male", "female", "male")),
    treatment = rep(c("A", "B"), each = 2)
  )

  out <- sumtab(~ age + sex | treatment, dat)

  expect_s3_class(out, "tbl_df")
  expect_equal(out$variable, c("age", "sex", "sex"))
  expect_equal(out$level[1], "")
  expect_true(all(c("Overall (n = 4)", "A (n = 2)", "B (n = 2)") %in% names(out)))
  expect_match(out[["Overall (n = 4)"]][1], "16.0")
  expect_match(out[["A (n = 2)"]][2], "1 \\(50.0%\\)")
})

test_that("sumtab adds one p-value and selected test per variable", {
  set.seed(1)
  dat <- tibble::tibble(
    y = c(stats::rnorm(40), stats::rnorm(40, 0.5)),
    x = rep(c("yes", "no"), 40),
    g = rep(c("A", "B"), each = 40)
  )

  out <- sumtab(~ y + x | g, dat, p_value = TRUE)

  expect_true(all(c("p.value", "test") %in% names(out)))
  expect_equal(out$test[out$variable == "y"][1], "Student independent t-test")
  expect_equal(out$test[out$variable == "x"][1], "Chi-square test of independence")
  expect_equal(out$p.value[out$variable == "x"][2], "")
  expect_equal(out$test[out$variable == "x"][2], "")
})

test_that("sumtab selects multi-group numeric tests", {
  set.seed(2)
  dat <- tibble::tibble(
    y = stats::rnorm(90),
    g = rep(c("A", "B", "C"), each = 30)
  )

  out <- sumtab(~ y | g, dat, p_value = TRUE)

  expect_equal(out$test[1], "One-way ANOVA")
  expect_false(is.na(out$p.value[1]))
})

test_that("sumtab selects Fisher exact test for sparse categorical tables", {
  dat <- tibble::tibble(
    x = c("yes", "yes", "yes", "no", "no", "yes"),
    g = c("A", "A", "B", "B", "B", "B")
  )

  out <- sumtab(~ x | g, dat, p_value = TRUE)

  expect_equal(out$test[1], "Fisher exact test")
  expect_false(is.na(out$p.value[1]))
})

test_that("sumtab validates formula shape and selected columns", {
  dat <- tibble::tibble(x = 1:3)

  expect_error(sumtab(x ~ 1, dat), "one-sided")
  expect_error(sumtab(~ missing, dat), "Missing required column")
})

test_that("sumtab prints without tibble type rows or missing-value markers", {
  dat <- tibble::tibble(
    age = c(10, 12, 20, 22),
    sex = factor(c("female", "male", "female", "male")),
    treatment = rep(c("A", "B"), each = 2)
  )

  txt <- capture.output(sumtab(~ age + sex | treatment, dat))

  expect_false(any(grepl("<chr>|<NA>", txt)))
  expect_true(any(grepl("Overall \\(n = 4\\)", txt)))
  expect_true(any(grepl("A \\(n = 2\\)", txt)))
})
