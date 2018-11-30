#' calcBatchEffectsForBatch
#'
#'@keywords internal
#'
#'@description  function to calculate batch effects for every gene in a batch
#'
#'@importFrom stats ks.test
#'@import BiocParallel
#'@import futile.logger
#'
#'@return the medians p-values for genes in a batch
calcBatchEffectsForBatch <- function(batch, samples, data, BPPARAM=bpparam()) {
    
    flog.debug(paste("Calculating the batch effect for batch", batch))
    
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
    
    flog.debug(paste("Calculating the p-values for batch", batch))
    p_values <- bpmapply(function(X,Y){
        if(all(is.na(X))){
            return(0.0)
        }else{
            return(ks.test(X, Y)$p.value)
        }}, X = features_batch, Y = features_other)
    
    flog.debug(paste("Calculating the median differences for batch", batch))
    medians <- bpmapply(function(X,Y){
        median(X, na.rm = TRUE) - median(Y, na.rm = TRUE)
    }, X = features_batch, Y = features_other)
    
    DF <- data.frame(p_values=unlist(p_values), medians=unlist(medians), 
                     row.names = unique(data$feature))
    
    return(DF)
}