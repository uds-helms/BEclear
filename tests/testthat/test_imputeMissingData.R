testthat::context("Testing the imputation of missing values")

testthat::test_that("No NAs, serial", {
  bpparam <- SerialParam(log = TRUE, threshold = "INFO")
  data <- as.matrix(Hilbert(50))

  set.seed(1, kind = "Mersenne-Twister")
  res <- imputeMissingData(
    data = data, BPPARAM = bpparam,
    rowBlockSize = 10, colBlockSize = 10, outputFormat = ""
  )
  colnames(res) <- NULL
  testthat::expect_equal(res, data)
})

testthat::test_that("NAs, serial", {
  bpparam <- SerialParam(log = TRUE, threshold = "INFO")
  data <- as.matrix(Hilbert(10))

  data_missing <- data
  data_missing[1, 3] <- NA
  data_missing[9, 5] <- NA

  set.seed(1, kind = "Mersenne-Twister")
  res <- imputeMissingData(
    data = data_missing, BPPARAM = bpparam,
    rowBlockSize = 0, colBlockSize = 0, outputFormat = ""
  )
  colnames(res) <- NULL
  data[1, 3] <- 0.116
  data[9, 5] <- 0.020
  testthat::expect_equal(res, data, tolerance = .01)
})
