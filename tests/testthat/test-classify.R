test_that("hess_classify classe les reactions", {
  expect_equal(hess_classify(-100), "Exothermique")
  expect_equal(hess_classify(100), "Endothermique")
  expect_equal(hess_classify(0), "Thermoneutre")
})
