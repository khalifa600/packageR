test_that("hess_reverse inverse une reaction", {
  r <- list(
    reactants = c(A = 1),
    products = c(B = 1),
    delta_h = -50
  )

  result <- hess_reverse(r)

  expect_equal(result$delta_h, 50)
  expect_equal(result$reactants, r$products)
  expect_equal(result$products, r$reactants)
})
