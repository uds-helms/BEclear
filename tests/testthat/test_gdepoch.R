testthat::context("Testing the gradient descent for a single epoch")

testthat::test_that("quadratic D matrix, only integer values given", {
  res1 <- gdepoch(
    L = matrix(c(1, 2)), R = t(matrix(c(3, 4))),
    D = matrix(c(1, 1, 1, 1), nrow = 2),
    lambda = 0.2, is = c(1, 2, 1, 2), js = c(1, 1, 2, 2),
    gamma = 0.01
  )

  testthat::expect_equal(res1$L, matrix(c(0.636, 1.132)), tolerance = .01)
  testthat::expect_equal(res1$R, t(matrix(c(2.748, 3.644))), tolerance = .01)
})

testthat::test_that("long D matrix, only integer values given", {
  res1 <- gdepoch(
    L = matrix(c(5, 6)), R = t(matrix(c(7, 4, 8))),
    D = matrix(c(7, 5, 3, 1, 3, 4), nrow = 2),
    lambda = 0.2, is = c(1, 2, 1, 2, 1, 2), js = c(1, 1, 2, 2, 3, 3),
    gamma = 0.01
  )

  testthat::expect_equal(res1$L, matrix(c(-6.22, -8.08)), tolerance = .01)
  testthat::expect_equal(res1$R, t(matrix(c(-0.268, -0.476, -1.01))), tolerance = .01)
})

testthat::test_that("wide D matrix, only integer values given", {
  res1 <- gdepoch(
    L = matrix(c(5, 1, 2)), R = t(matrix(c(3, 3))),
    D = matrix(c(5, 6, 10, 11, 3, 1), nrow = 3),
    lambda = 0.2, is = c(1, 2, 3, 1, 2, 3), js = c(1, 1, 1, 2, 2, 2),
    gamma = 0.01
  )

  testthat::expect_equal(res1$L, matrix(c(4.14, 1.176, 1.93)), tolerance = .01)
  testthat::expect_equal(res1$R, t(matrix(c(2.21, 2.39))), tolerance = .01)
})
