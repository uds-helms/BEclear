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
    
    # flog.debug(paste("Calculating the medians for batch", batch))
    # 
    # median_batch <- samples[batch_id == batch][data, , 
    #                                            on=.(sample_id = sample), nomatch=0][, 
    #                                                                                 median(beta.value, 
    #                                                                                        na.rm = TRUE),
    #                                                                                 by=.(feature)]$V1
    # median_others <- samples[batch_id != batch][data, , 
    #                                             on=.(sample_id = sample), nomatch=0][, 
    #                                                                                  median(beta.value, 
    #                                                                                         na.rm = TRUE),
    #                                                                                  by=.(feature)]$V1
    # DF <- data.frame(median_batch - median_others, 
    #                  row.names = rownames(median_batch))
    # colnames(DF) <- batch
    # return(DF)
    
    flog.debug(paste("Calculating the medians for batch", batch))

    flog.debug("Subsetting data.table for batch")
    DT_batch <- samples[batch_id == batch][data,
                                           , on=.(sample_id = sample),
                                           nomatch=0][, .(feature, beta.value)]
    flog.debug("Subsetting data.table for others")
    DT_other <- samples[batch_id != batch][data,
                                           , on=.(sample_id = sample),
                                           nomatch=0][, .(feature, beta.value)]

    features <- lapply(unique(data$feature), FUN=function(x, DT_batch, DT_other)
    {
        return(list(batch=DT_batch[feature==x, beta.value],
                    others=DT_other[feature==x, beta.value]))
    }, DT_batch = DT_batch, DT_other = DT_other)



    medians <- bplapply(features, FUN = function(x)
        {
        return( median(x$batch, na.rm = TRUE) - median(x$others, na.rm = TRUE))
        }, BPPARAM=BPPARAM)

    DF <- data.frame(unlist(medians), row.names = unique(data$feature))
    colnames(DF) <- batch
    return(DF)
        
    }