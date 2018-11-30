testthat::context("Testing the detection and correction of Batch Effects")

testthat::test_that("3 batches, NAs, result including values below 0",{
    
    data <- data.table(sample=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5), 
                       feature=c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3),
                       beta.value=c(NA,  0.08099936,  0.18366184,  
                                    0.16625504, -1.26959907,  NA, 
                                    -1.41200541, -0.01696149, -0.54431935,
                                    -0.6548353,  NA, -1.2723697,
                                    0.5380646,  NA, -1.1657070))
    
    samples<-data.table(sample_id = as.character(c(1,2,3,4,5)), batch_id=c(1,1,2,3,3))
    DT<-as.matrix(dcast(data = data, formula = feature ~ sample, 
                        value.var = "beta.value"))[,2:6]
    row.names(DT)<-c("1", "2", "3")
    
    DT_expected <- matrix(c(0.02496767,0.08099936, 0.18366184,
                            0.16625504, 0, 0.02423866,
                            0, 0, 0,
                            0, 0.01631574, 0,
                            0.5380646,0, 0),ncol = 5)
    colnames(DT_expected) <- c("1", "2", "3", "4", "5")
    row.names(DT_expected) <-c("1", "2", "3")
    
    res<-correctBatchEffect(data=DT, samples = samples, outputFormat = "")
    testthat::expect_equal(res$correctedPredictedData, DT_expected, tolerance = .001)
}
)