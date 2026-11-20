test_that("hess_combine combine les reactions", {
  r1 <- list(
    reactants = c(A = 1),
    products = c(B = 1),
    delta_h = -50
  )

  r2 <- list(
    reactants = c(B = 1),
    products = c(C = 1),
    delta_h = -20
  )

  result <- hess_combine(list(r1, r2))

  expect_equal(result$delta_h, -70)
  expect_equal(result$reactants, c(A = 1))
  expect_equal(result$products, c(C = 1))
})

test_that("hess_combine accepte des coefficients negatifs", {
  r1 <- list(reactants = c(A = 1), products = c(B = 1), delta_h = -50)
  result <- hess_combine(list(r1), coefficients = -1)

  expect_equal(result$delta_h, 50)
  expect_equal(result$reactants, c(B = 1))
  expect_equal(result$products, c(A = 1))
})
