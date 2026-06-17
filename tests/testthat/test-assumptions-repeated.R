test_that("repeated workflow reports sphericity note", {
  dat <- make_assumption_data(60)
  x <- test_repeated(dat, c(sbp_baseline, sbp_3m, sbp_6m), id = id, plot = FALSE)
  expect_true(any(grepl("Sphericity", x$assumptions$name)))
})
