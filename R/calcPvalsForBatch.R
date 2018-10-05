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
    
    flog.debug("Performing a Kolmogorov-Smirnov tests to obtain the p-value")
    
    p_values <- vapply(X=unique(data$feature), FUN=function(f, x, y){
        x_bat <- x[feature == f, beta.value]
        if(all(is.na(x_bat))){
            0.0
        }else{
            ks.test(x_bat, y[feature == f, beta.value])$p.value
        }
        
    }, FUN.VALUE = numeric(1), x = DT_batch, y = DT_other)

   
    DF <- data.frame(p_values, row.names = unique(data$feature))
    colnames(DF) <- batch
    return(DF)
}