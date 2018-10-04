#' calcPvalsForBatch
#'
#'@keywords internal
#'
#'@description  function to calculate p-values for every gene in a batch
#'
#'@importFrom stats ks.test
#'
#'@return the p-values for genes in a batch
calcPvalsForBatch <- function(batch, samples, data) {
    
    DT_batch <- samples[batch_id == batch][data,
                  , on=.(sample_id = sample), 
                  nomatch=0][, .(feature, beta.value)]
    DT_other <- samples[batch_id == batch][data,
                                           , on=.(sample_id = sample), 
                                           nomatch=0][, .(feature, beta.value)]
    
    suppressWarnings(
        p_values <- DT_batch[ ,if(all(is.na(beta.value))){0.0} 
                              else{ks.test(beta.value, 
                                           DT_other[feature==feature, 
                                                    beta.value], 
                                           exact = FALSE)$p.value}, 
                              by=.(feature)])
    DF <- data.frame(p_values$V1, row.names = p_values$feature)
    colnames(DF) <- batch
    return(DF)
}