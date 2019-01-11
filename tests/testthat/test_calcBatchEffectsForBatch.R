testthat::context("Testing the caculations of batch effects for single batch")

testthat::test_that("2 batches, only available values", {
  data <- data.table(
    sample = c(1, 1, 1, 2, 2, 2, 3, 3, 3), feature = c(1, 2, 3, 1, 2, 3, 1, 2, 3),
    beta.value = c(
      0.29382666, 0.08099936, 0.18366184,
      0.16625504, -1.26959907, 2.34949332,
      -1.41200541, -0.01696149, -0.54431935
    )
  )

  samples <- data.table(sample_id = c(1, 2, 3), batch_id = c(1, 1, 2))

  res1 <- calcBatchEffectsForBatch(batch = 1, samples = samples, data = data,
                                   BPPARAM = SerialParam())
  res2 <- calcBatchEffectsForBatch(batch = 2, samples = samples, data = data,
                                   BPPARAM = SerialParam())

  testthat::expect_equal(res1[, 1], c(0.667, 1.000, 0.667), tolerance = .001)
  testthat::expect_equal(res2[, 1], c(0.667, 1.000, 0.667), tolerance = .001)
  testthat::expect_equal(res1[, 2], c(1.642, -0.577, 1.811), tolerance = .001)
  testthat::expect_equal(res2[, 2], c(-1.642, 0.577, -1.811), tolerance = .001)
})

testthat::test_that("3 batches, NAs", {
  data <- data.table(
    sample = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5),
    feature = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3),
    beta.value = c(
      NA, 0.08099936, 0.18366184,
      0.16625504, -1.26959907, NA,
      -1.41200541, -0.01696149, -0.54431935,
      -0.6548353, NA, -1.2723697,
      0.5380646, NA, -1.1657070
    )
  )

  samples <- data.table(sample_id = c(1, 2, 3, 4, 5), batch_id = c(1, 1, 2, 3, 3))

  res1 <- calcBatchEffectsForBatch(batch = 1, samples = samples, data = data,
                                   BPPARAM = SerialParam())
  res2 <- calcBatchEffectsForBatch(batch = 2, samples = samples, data = data,
                                   BPPARAM = SerialParam())
  res3 <- calcBatchEffectsForBatch(batch = 3, samples = samples, data = data,
                                   BPPARAM = SerialParam())

  testthat::expect_equal(res1[, 1], c(1.00, 1.00, 0.50), tolerance = .001)
  testthat::expect_equal(res2[, 1], c(0.50, 1.00, 1.00), tolerance = .001)
  testthat::expect_equal(res3[, 1], c(1.00, 0.00, 0.333), tolerance = .001)

  testthat::expect_equal(res1[, 2], c(0.821, -0.577, 1.349), tolerance = .001)
  testthat::expect_equal(res2[, 2], c(-1.578, 0.577, 0.621), tolerance = .001)
  testthat::expect_equal(res3[, 2], c(0.564, NA, -1.039), tolerance = .001)
})
