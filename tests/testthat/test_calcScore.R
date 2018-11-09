testthat::context("Testing the calculation of BE scores")

testthat::test_that("matrix with batch effects",{
    data(BEclearData)
    ex.data <- ex.data[31:90,7:26]
    ex.samples <- ex.samples[7:26,]
    library(data.table)
    samples <- data.table(ex.samples)
    data <- data.table(feature=rownames(ex.data), ex.data)
    data <- melt(data = data, id.vars = "feature", variable.name = "sample",
                value.name = "beta.value")
    setkey(data, "feature", "sample")
    med <- calcMedians(data=data, samples=samples)
    pvals <- calcPvalues(data=data, samples=samples, adjusted=TRUE,
    method="fdr")
    sum <- calcSummary(medians=med, pvalues=pvals)
    score.table <- calcScore(data=ex.data, samples=ex.samples, summary=sum,
                             saveAsFile=FALSE)
    
    testthat::expect_equal(score.table$BEscore, c(0,0,0,0.73333333))
    
}
)