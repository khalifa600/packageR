test_that("hess_multiply multiplie correctement", {
  r <- list(
    reactants = c(A = 1),
    products = c(B = 1),
    delta_h = -50
  )

  result <- hess_multiply(r, 2)

  expect_equal(result$delta_h, -100)
  expect_equal(result$reactants, c(A = 2))
})

test_that("hess_multiply refuse un facteur negatif", {
  r <- list(reactants = c(A = 1), products = c(B = 1), delta_h = -50)
  expect_error(hess_multiply(r, -1))
})
