#' calcPvalsForBatch
#'
#'@keywords internal
#'
#'@description  function to calculate p-values for every gene in a batch
#'
#'@importFrom stats ks.test
#'@import futile.logger
#'
#'@return the p-values for genes in a batch
calcPvalsForBatch <- function(batch, samples, data) {
    
    flog.debug(paste("Calculating the pvalues for batch", batch))
    
    flog.debug("Subsetting data.table for batch")
    DT_batch <- samples[batch_id == batch][data,
                  , on=.(sample_id = sample), 
                  nomatch=0][, .(feature, beta.value)]
    flog.debug("Subsetting data.table for others")
    DT_other <- samples[batch_id != batch][data,
                                           , on=.(sample_id = sample), 
                                           nomatch=0][, .(feature, beta.value)]
    
    flog.debug("Performing a Kolmogorov-Smirnov test to obtain the p-value")
    
    p_values <- vector(mode = "numeric", length = length(unique(data$feature)))
    i <- 1
    for (f in unique(data$feature)){
        flog.debug(paste("Calculating for feature", f))
        batch_betas <- DT_batch[feature == f, beta.value]
        if(all(is.na(batch_betas))){
            p_values[i] <- 0.0
        }else{
            p_values[i] <- ks.test(batch_betas,
                                   DT_other[feature == f, beta.value])$p.value
        }
        
        
        i <- i + 1
    }
    
    # suppressWarnings(
    #     p_values <- DT_batch[ ,if(all(is.na(beta.value))){0.0} 
    #                           else{ks.test(beta.value, 
    #                                        DT_other[feature==feature, 
    #                                                 beta.value], 
    #                                        exact = FALSE)$p.value}, 
    #                           by=.(feature)])
    #DF <- data.frame(p_values$V1, row.names = p_values$feature)
    DF <- data.frame(p_values, row.names = unique(data$feature))
    colnames(DF) <- batch
    return(DF)
}