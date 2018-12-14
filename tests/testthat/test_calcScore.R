testthat::context("Testing the calculation of BE scores")

testthat::test_that("matrix with batch effects", {
  data(BEclearData)
  ex.data <- ex.data[31:90, 7:26]
  ex.samples <- ex.samples[7:26, ]
  batchEffects <- calcBatchEffects(data = ex.data, samples = ex.samples,
                                    adjusted = TRUE, method = "fdr")
  med <- batchEffects$med
  pvals <- batchEffects$pval
  
  
  sum <- calcSummary(medians = med, pvalues = pvals)
  score.table <- calcScore(
    data = ex.data, samples = ex.samples, summary = sum,
    saveAsFile = FALSE
  )

  testthat::expect_equal(score.table$BEscore, c(0, 0, 0, 0.73333333))
})
