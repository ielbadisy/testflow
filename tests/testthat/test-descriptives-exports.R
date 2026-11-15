test_that("descriptivesnumeric matches base summaries", {
  d <- descriptivesnumeric(mtcars, "mpg", group = "am")
  expect_equal(nrow(d), 2)
  expect_equal(sort(unname(d$mean)), sort(as.vector(tapply(mtcars$mpg, mtcars$am, mean))))
  expect_error(descriptivesnumeric(mtcars, "nope"))
})

test_that("descriptivescategorical matches table()", {
  d <- descriptivescategorical(mtcars, "cyl")
  expect_equal(sum(d$n), nrow(mtcars))
  expect_equal(sort(unname(d$n)), sort(as.vector(table(mtcars$cyl))))
})
