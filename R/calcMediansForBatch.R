#' calcMediansForBatch
#'
#' @keywords internal
#'
#' @description  function to calculate medians for every gene for a batch
#' 
#' @importFrom stats median
#' @import BiocParallel
#' @import data.table
#' @import futile.logger
#'
#' @return the medians for genes in a batch
calcMediansForBatch <- function(batch, samples, data, BPPARAM=bpparam()) {
    

    
    flog.debug(paste("Calculating the medians for batch", batch))

    flog.debug("Subsetting data.table for batch")
    DT_batch <- samples[batch_id == batch][data,
                                           , on=.(sample_id = sample),
                                           nomatch=0][, .(feature, beta.value)]
    flog.debug("Subsetting data.table for others")
    DT_other <- samples[batch_id != batch][data,
                                           , on=.(sample_id = sample),
                                           nomatch=0][, .(feature, beta.value)]

    features_batch <- DT_batch[, list(list(beta.value)), by=feature]$V1
    features_other <- DT_other[, list(list(beta.value)), by=feature]$V1
    
    features<-mapply(function(X,Y){list(batch = X, others = Y)}, 
                     X = features_batch, Y = features_other, SIMPLIFY = F)



    medians <- bplapply(features, FUN = function(x)
        {
        return( median(x$batch, na.rm = TRUE) - median(x$others, na.rm = TRUE))
        }, BPPARAM=BPPARAM)

    DF <- data.frame(unlist(medians), row.names = unique(data$feature))
    colnames(DF) <- batch
    return(DF)
        
    }