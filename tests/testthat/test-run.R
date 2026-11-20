test_that("hess_run execute une analyse", {
  r <- list(
    reactants = c(A = 1),
    products = c(B = 1),
    delta_h = -100
  )

  result <- hess_run(
    reactions = list(r),
    show_diagram = FALSE,
    verbose = FALSE
  )

  expect_equal(result$delta_h, -100)
  expect_equal(result$classification, "Exothermique")
})

test_that("hess_run trouve seul les coefficients a partir d'une cible", {
  ex <- hess_example("formation_methane")

  result <- hess_run(
    ex$reactions,
    target = ex$target,
    show_diagram = FALSE,
    verbose = FALSE
  )

  expect_equal(result$delta_h, -74.8)
  expect_equal(as.numeric(result$coefficients), c(1, 2, -1))
  expect_true(result$target_reached)
  expect_true(result$atom_check$balanced)
  expect_equal(result$reactants, c("C(s)" = 1, "H2(g)" = 2))
  expect_equal(result$products, c("CH4(g)" = 1))
})

test_that("hess_run fonctionne pour le monoxyde de carbone", {
  ex <- hess_example("monoxyde_carbone")
  result <- hess_run(ex$reactions, target = ex$target,
                     show_diagram = FALSE, verbose = FALSE)

  expect_equal(result$delta_h, -110.5)
})

test_that("hess_run accepte une reaction seule et produit un diagramme", {
  skip_if_not_installed("ggplot2")

  ex <- hess_example("methane")
  pdf(NULL)
  on.exit(dev.off())
  result <- hess_run(ex$reactions[[1]], verbose = FALSE)

  expect_equal(result$delta_h, -890.3)
  expect_s3_class(result$diagram, "ggplot")
})

test_that("hess_run signale une cible impossible", {
  r <- list(reactants = c(A = 1), products = c(B = 1), delta_h = -10)
  target <- list(reactants = c(C = 1), products = c(D = 1))

  expect_error(
    hess_run(list(r), target = target, show_diagram = FALSE, verbose = FALSE)
  )
})
