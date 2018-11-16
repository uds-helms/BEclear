#' calcPvalsForBatch
#'
#'@keywords internal
#'
#'@description  function to calculate p-values for every gene in a batch
#'
#'@importFrom stats ks.test
#'@import BiocParallel
#'@import futile.logger
#'
#'@return the p-values for genes in a batch
calcPvalsForBatch <- function(batch, samples, data, BPPARAM=bpparam()) {
    
    flog.debug(paste("Calculating the pvalues for batch", batch))
    
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
    
    
    
    p_values <- bplapply(features, calcPvalsForGene, BPPARAM=BPPARAM)

    
    DF <- data.frame(unlist(p_values), row.names = unique(data$feature))
    colnames(DF) <- batch
    return(DF)
}