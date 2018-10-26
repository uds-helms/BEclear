testthat::context("Replacement of values below 0 or above 1")

testthat::test_that("Matrix with values above and below",{
    res1 <- replaceWrongValues(matrix(c(2, 0.5, 0.4, 0.2, -1, -0.1), ncol = 2))
    
    testthat::expect_equal(as.numeric(res1[,1]), c(1.0, 0.5, 0.4))
    testthat::expect_equal(as.numeric(res1[,2]), c(0.2, 0.0, 0.0))
}
)