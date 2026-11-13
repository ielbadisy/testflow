test_that("paired workflow checks paired differences", {
  dat <- .bt(before = rnorm(40), after = rnorm(40))
  x <- testpaired(dat, after ~ before, plot = FALSE)
  expect_true(any(grepl("paired differences", x$assumptions$name, ignore.case = TRUE)))
  expect_false(any(grepl("^Normality by group$", x$assumptions$name)))
})
