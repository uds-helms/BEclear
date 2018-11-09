testthat::context("Testing the summary calculation")

testthat::test_that("matrix with batch effects",{
    sum <- calcSummary(data.frame( b1=c(0.9, NA, 0.4) , b2=c(NA, 0.001, 0.8), 
                                   b3=c(0.2, 0.5, 0.01)), 
                       data.frame( b1=c(0.9, 0.3, 0.04) , b2=c(0.01, 0.002, 0.06), 
                                   b3=c(0.9, 0.00001, 0.4)))
    testthat::expect_equal(sum$median, c(NA, 0.5))
    testthat::expect_equal(sum$pvalue, c(0.01, 0.00001))
    
}
)

testthat::test_that("matrix without batch effects",{
    sum <- calcSummary(data.frame( b1=c(0.1, 0.2, 0.4) , b2=c(0.01, 0.001, 0.8), 
                                   b3=c(0.2, 0.5, 0.01)), 
                       data.frame( b1=c(0.9, 0.3, 0.04) , b2=c(0.01, 0.002, 0.06), 
                                   b3=c(0.9, 0.1, 0.4)))
    
    testthat::expect_null(sum)
}
)