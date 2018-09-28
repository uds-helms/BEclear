#' calcMediansForBatch
#'
#' @keywords internal
#'
#' @description  function to calculate medians for every gene for a batch
#' 
#' @importFrom stats median
#' @import data.table
#' @import futile.logger
#'
#' @return the medians for genes in a batch
calcMediansForBatch <- function(batch, samples, data) {
    
    flog.debug(paste("Calculating the medians for batch", batch))
    
    median_batch <- samples[batch_id == batch][data, , 
                         on=.(sample_id = sample), nomatch=0][, 
                                                   median(beta.value, 
                                                          na.rm = TRUE),
                                                   by=.(feature)]$V1
    median_others <- samples[batch_id != batch][data, , 
                          on=.(sample_id = sample), nomatch=0][, 
                                                    median(beta.value, 
                                                           na.rm = TRUE),
                                                    by=.(feature)]$V1
    DF <- data.frame(median_batch - median_others, 
                     row.names = rownames(median_batch))
    colnames(DF) <- batch
    return(DF)
        
    }