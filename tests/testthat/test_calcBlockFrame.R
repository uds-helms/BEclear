testthat::context("Calculation of block positions")

testthat::test_that("Values that add up",{
    res1<-calcBlockFrame(c(2,1,11,10,20),c(2,1,11,10,20),10,10)
    
    testthat::expect_equal(res1[,1], c(1,2,3,4))
    testthat::expect_equal(res1[,2], c(1,1,11,11))
    testthat::expect_equal(res1[,3], c(10,10,20,20))
    testthat::expect_equal(res1[,4], c(1,11,1,11))
    testthat::expect_equal(res1[,5], c(10,20,10,20))
}
)

testthat::test_that("Values that don't add up",{
    res1<-calcBlockFrame(c(2,1,11,10,19),c(2,1,11,10,19),11,11)
    
    testthat::expect_equal(res1[,1], c(1,2,3,4))
    testthat::expect_equal(res1[,2], c(1,1,11,11))
    testthat::expect_equal(res1[,3], c(10,10,19,19))
    testthat::expect_equal(res1[,4], c(1,11,1,11))
    testthat::expect_equal(res1[,5], c(10,19,10,19))
}
)