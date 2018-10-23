testthat::context("Counting of non zero revealed entries")

testthat::test_that("Only existing indices",{
    testthat::expect_equal( nnz(c(3,4,5,2,1), 5), c(1,1,1,1,1))
    testthat::expect_equal( nnz(c(1,3,1,2,1), 3), c(3,1,1))
    testthat::expect_equal( nnz(c(2,2,1,2,1), 2), c(2,3))
}
)

testthat::test_that("existing and missing indices",{
    testthat::expect_equal( nnz(c(3,4,5,2,1), 6), c(1,1,1,1,1, 0))
    testthat::expect_equal( nnz(c(2,2,2,2,2), 3), c(0, 5, 0))
}
)