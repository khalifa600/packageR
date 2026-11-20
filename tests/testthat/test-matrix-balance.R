test_that("hess_matrix construit la matrice stoechiometrique", {
  ex <- hess_example("formation_methane")
  m <- hess_matrix(ex$reactions, ex$target)

  expect_equal(ncol(m$matrix), 3)
  expect_equal(nrow(m$matrix), length(m$species))
  expect_equal(length(m$target), length(m$species))
  expect_equal(m$matrix["CH4(g)", 3], -1)
  expect_equal(m$target[["CH4(g)"]], 1)
})

test_that("hess_check_balance verifie les atomes", {
  ex <- hess_example("methane")
  expect_true(hess_check_balance(ex$reactions[[1]])$balanced)

  bad <- list(
    reactants = c("CH4(g)" = 1, "O2(g)" = 1),
    products = c("CO2(g)" = 1, "H2O(l)" = 2),
    delta_h = -1
  )
  expect_false(hess_check_balance(bad)$balanced)

  ca <- list(
    reactants = c("Ca(OH)2(s)" = 1),
    products = c("CaO(s)" = 1, "H2O(l)" = 1),
    delta_h = 65
  )
  expect_true(hess_check_balance(ca)$balanced)
})
