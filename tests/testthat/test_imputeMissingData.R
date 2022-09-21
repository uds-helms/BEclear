testthat::context("Testing the imputation of missing values")


testthat::test_that("No NAs, serial", {
  library(seewave)
  bpparam <- SerialParam(log = TRUE, threshold = "INFO")
  data <- as.matrix(Hilbert(50))

  set.seed(1, kind = "Mersenne-Twister")
  res <- imputeMissingData(
    data = data, BPPARAM = bpparam,
    rowBlockSize = 10, colBlockSize = 10, outputFormat = ""
  )
  colnames(res) <- NULL
  testthat::expect_equal(object = as.vector(res), expected =  as.vector(data))
})

testthat::test_that("NAs, serial", {
  library(seewave)
  bpparam <- SerialParam(log = TRUE, threshold = "INFO")
  data <- as.matrix(Hilbert(50))

  data_missing <- data
  data_missing[1, 3] <- NA
  data_missing[9, 5] <- NA

  set.seed(1, kind = "Mersenne-Twister")
  res <- imputeMissingData(
    data = data_missing, BPPARAM = bpparam,
    rowBlockSize = 10, colBlockSize = 10, outputFormat = ""
  )
  colnames(res) <- NULL
  data[1, 3] <- 0.116
  data[9, 5] <- 0.020
  testthat::expect_equal(object = as.vector(res), expected =  as.vector(data), 
                         tolerance = .01)
})
