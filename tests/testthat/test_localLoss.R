testthat::context("Testing the local loss")

testthat::test_that("quadratic error matrix, only integer values given",{
    
    res1<-localLoss(L = matrix(c(1,2)), R = t(matrix(c(3,4))), 
                  error_matrix = matrix(c(1,1,1,1), nrow = 2),
                  is = c(1,2,1,2), js = c(1,1,2,2))
    
    testthat::expect_equal(res1$dL, matrix(c(7, 7)), tolerance = .01)
    testthat::expect_equal(res1$dR, t(matrix(c(3, 3))), tolerance = .01)
}
)

testthat::test_that("long error matrix, only integer values given",{
    
    res1<-localLoss(L = matrix(c(5,6)), R = t(matrix(c(7,4,8))), 
                  error_matrix = matrix(c(7,5,3,1,3,4), nrow = 2),
                  is = c(1,2,1,2,1,2), js = c(1,1,2,2,3,3))
    
    testthat::expect_equal(res1$dL, matrix(c(85, 71)), tolerance = .01)
    testthat::expect_equal(res1$dR, t(matrix(c(65, 21, 39))), tolerance = .01)
}
)

testthat::test_that("wide error matrix, only integer values given",{
    
    res1<-localLoss(L = matrix(c(5,1,2)), R = t(matrix(c(3,3))), 
                  error_matrix = matrix(c(5,6,10,11,3,1), nrow = 3),
                  is = c(1,2,3,1,2,3), js = c(1,1,1,2,2,2))
    
    testthat::expect_equal(res1$dL, matrix(c(48, 27,33)), tolerance = .01)
    testthat::expect_equal(res1$dR, t(matrix(c(51, 60))), tolerance = .01)
}
)