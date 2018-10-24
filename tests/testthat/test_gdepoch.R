testthat::context("Testing the gradient descent for a single epoch")

testthat::test_that("matrices with integer values given",{
    
    res1<-gdepoch(L = matrix(c(1,2)), R = t(matrix(c(3,4))), 
                  D = matrix(c(1,1,1,1), nrow = 2),
                  lambda = 0.2, is = c(1,2,1,2), js = c(1,1,2,2), 
                  nnzis = c(2,2), nnzjs = c(2,2), gamma = 0.01)
    
    testthat::expect_equal(res1$L, matrix(c(0.636, 1.132)), tolerance = .001)
    testthat::expect_equal(res1$R, t(matrix(c(2.748, 3.644))), tolerance = .001)
}
)