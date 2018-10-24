testthat::context("calculating the loss for matrices L and R")

testthat::test_that("matrices with integer results",{
    res1<-loss(L = matrix(c(1,2)), R = t(matrix(c(3,4))), 
         D = matrix(c(1,1,1,1), nrow = 2),
         lambda = 0.2)
    testthat::expect_equal( res1$loss, 93)
    testthat::expect_equal( res1$error_matrix, matrix(c(-2, -5, -3, -7), ncol = 2))
}
)