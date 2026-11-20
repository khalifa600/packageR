test_that("hess_calculate calcule correctement", {
  expect_equal(
    hess_calculate(c(1, 2), c(-100, -50)),
    -200
  )
})

test_that("hess_calculate refuse les entrees invalides", {
  expect_error(hess_calculate(c(1, 2), c(-100)))
  expect_error(hess_calculate("a", 1))
  expect_error(hess_calculate(NA_real_, 1))
})
