testthat::context("Testing the detection and correction of Batch Effects")

testthat::test_that("3 batches, NAs", {
  data <- data.table(
    sample = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    feature = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3),
    beta.value = c(
      NA, 0.08099936, 0.18366184,
      0.16625504, 0.26959907, NA,
      0.41200541, 0.01696149, 0.54431935,
      0.6548353, NA, 0.2723697,
      0.5380646, NA, 0.1657070
    )
  )

  samples <- data.table(sample_id = as.character(c(1, 2, 3, 4, 5)), batch_id = c(1, 1, 2, 3, 3))
  DT <- as.matrix(dcast(
    data = data, formula = feature ~ sample,
    value.var = "beta.value"
  ))[, 2:6]
  row.names(DT) <- c("1", "2", "3")

  DT_expected <- matrix(c(
    0.01401, 0.08099936, 0.18366184,
    0.16625504, 0.26960, 0.01138,
    0.41201, 0.01696, 0.54432,
    0.65484, 0.00615, 0.27237,
    0.5380646, 0.0048, 0.1657
  ), ncol = 5)
  colnames(DT_expected) <- c("1", "2", "3", "4", "5")
  row.names(DT_expected) <- c("1", "2", "3")

  set.seed(1, kind = "Mersenne-Twister")
  res <- correctBatchEffect(data = DT, samples = samples, outputFormat = "")
  testthat::expect_equal(res$correctedPredictedData, DT_expected, tolerance = .001)
})

# testthat::test_that("3 batches, invalid values", {
#   samples <- data.table(sample_id = as.character(c(1, 2, 3, 4, 5)), batch_id = c(1, 1, 2, 3, 3))
#   DT <- matrix(c(
#     0.01401, 0.08099936, 0.18366184,
#     2, 0.26960, 0.01138,
#     0.41201, 0.01696, -0.54432,
#     0.65484, 0.00615, 0.27237,
#     0.5380646, 6, 0.1657
#   ), ncol = 5)
#   row.names(DT) <- c("1", "2", "3")
#   colnames(DT) <- c("1", "2", "3", "4", "5")
# 
#   DT_expected <- matrix(c(
#     0.01401, 0.08099936, 0.18366184,
#     0.00163, 0.26960, 0.01138,
#     0.41201, 0.01696, 0.01079,
#     0.65484, 0.00615, 0.27237,
#     0.5380646, 0.00167, 0.1657
#   ), ncol = 5)
#   colnames(DT_expected) <- c("1", "2", "3", "4", "5")
#   row.names(DT_expected) <- c("1", "2", "3")
# 
#   set.seed(1, kind = "Mersenne-Twister")
#   res <- correctBatchEffect(data = DT, samples = samples, outputFormat = "")
#   testthat::expect_equal(res$correctedPredictedData, DT_expected, tolerance = .001)
# })

testthat::test_that("3 batches, duplicated colnames ", {
  samples <- data.table(sample_id = as.character(c(1, 2, 3, 4, 5)), batch_id = c(1, 1, 2, 3, 3))
  DT <- matrix(c(
    0.91401, 0.8099936, 0.18366184,
    2, 0.26960, 0.001138,
    0.41201, 0.1696, 0.54432,
    0.65484, 0.00615, 0.27237,
    0.5380646, 6, 0
  ), ncol = 5)
  row.names(DT) <- c("1", "2", "3")
  colnames(DT) <- c("1", "1", "1", "4", "5")

  DT_expected <- matrix(c(
    0.91401, 0.80999, 0.18366184,
    0.05021, 0.26960, 0.00114,
    0.41201, 0.16960, 0.54432,
    0.65484, 0.00615, 0.27237,
    0.5380646, 0.08527, 0
  ), ncol = 5)
  colnames(DT_expected) <- c("1", "2", "3", "4", "5")
  row.names(DT_expected) <- c("1", "2", "3")

  set.seed(1, kind = "Mersenne-Twister")
  res <- correctBatchEffect(data = DT, samples = samples, outputFormat = "")
  testthat::expect_equal(res$correctedPredictedData, DT_expected, tolerance = .001)
})
