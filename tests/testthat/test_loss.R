testthat::context("calculating the loss for matrices L and R")

testthat::test_that("matrices with integer results",{
    testthat::expect_equal( loss(L = matrix(c(1,2)), R = t(matrix(c(3,4))), 
                            D = matrix(c(1,1,1,1), nrow = 2),
                            lambda = 0.2, is = c(1,2,1,2), 
                            js = c(1,1,2,2)), 93)
}
)