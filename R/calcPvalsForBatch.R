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
    .Deprecated("calcBatchEffects")
    
    flog.debug(paste("Calculating the pvalues for batch", batch))
    
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
    
    
    p_values <- bpmapply(function(X,Y){
        if(all(is.na(X))){
            return(0.0)
            }else{
            return(ks.test(X, Y)$p.value)
        }}, X = features_batch, Y = features_other)

    
    DF <- data.frame(unlist(p_values), row.names = unique(data$feature))
    colnames(DF) <- batch
    return(DF)
}