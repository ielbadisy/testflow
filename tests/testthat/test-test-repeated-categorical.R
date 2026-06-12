test_that("test_repeated_categorical computes true Cochran Q", {
  mat <- matrix(
    c(
      1, 0, 1,
      1, 1, 1,
      0, 0, 1,
      1, 0, 0,
      0, 1, 1,
      1, 1, 0
    ),
    ncol = 3,
    byrow = TRUE
  )
  dat <- tibble::tibble(t1 = mat[, 1], t2 = mat[, 2], t3 = mat[, 3])
  x <- test_repeated_categorical(dat, c(t1, t2, t3))
  col_totals <- colSums(mat)
  row_totals <- rowSums(mat)
  k <- ncol(mat)
  q <- ((k - 1) * (k * sum(col_totals^2) - sum(col_totals)^2)) /
    (k * sum(row_totals) - sum(row_totals^2))
  expect_s3_class(x, "testflow_repeated_categorical")
  expect_equal(x$recommended$test, "Cochran Q test")
  expect_equal(unname(x$primary_test$statistic[1]), q)
  expect_equal(x$primary_test$p.value[1], stats::pchisq(q, df = 2, lower.tail = FALSE))
  expect_true(inherits(x$posthoc, "data.frame"))
  expect_equal(unique(x$posthoc$method), "McNemar test")
})
