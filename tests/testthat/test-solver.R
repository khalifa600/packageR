test_that("hess_solver trouve les coefficients", {
  matrice <- matrix(
    c(1, 0, 0, 1),
    nrow = 2
  )

  result <- hess_solver(
    matrice,
    c(2, 3)
  )

  expect_equal(
    as.numeric(result$coefficients),
    c(2, 3)
  )

  expect_true(result$exact)
})

test_that("hess_solver gere un systeme sous-determine", {
  matrice <- matrix(c(1, 1), nrow = 1)
  result <- hess_solver(matrice, 2)

  expect_true(result$exact)
  expect_false(result$unique)
  expect_equal(sum(result$coefficients), 2)
})

test_that("hess_solver signale une cible impossible", {
  matrice <- matrix(c(1, 1), nrow = 2)
  expect_warning(result <- hess_solver(matrice, c(1, 2)))
  expect_false(result$exact)
})
