testthat::context("Calculation of block positions")

testthat::test_that("Values that add up",{
    testthat::expect_equal(calcPositions(100,10), c(10,81,91,90,100))
    testthat::expect_equal(calcPositions(27,3), c(9,22,25,24,27))
    testthat::expect_equal(calcPositions(16,4), c(4,9,13,12,16))
}
)


testthat::test_that("Values that don't add up",{
    testthat::expect_equal(calcPositions(99,10), c(10,81,91,90,99))
    testthat::expect_equal(calcPositions(200,13), c(15,170,183,182,200))
    testthat::expect_equal(calcPositions(1363,13), c(105,1340,1353,1352,1363))
}
)

testthat::test_that("Block size smaller than 1", {
    testthat::expect_error(calcPositions(5,-1))
    testthat::expect_error(calcPositions(33,0))
    testthat::expect_error(calcPositions(1,-3))
}
)

testthat::test_that("Number of entries smaller than 1", {
    testthat::expect_error(calcPositions(0,1))
    testthat::expect_error(calcPositions(0,0))
    testthat::expect_error(calcPositions(-1,-3))
    testthat::expect_error(calcPositions(-20,67))
}
)

testthat::test_that("Number smallerthan block size", {
    testthat::expect_equal(calcPositions(13,55), c(1,0,1,0,13))

}
)