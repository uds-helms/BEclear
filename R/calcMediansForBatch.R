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
    
    median_batch <- data[samples, , 
                         on=.(sample = sample_id)][batch_id == batch, 
                                                   median(beta.value, 
                                                          na.rm = TRUE),
                                                   by=.(feature)]$V1
    median_others <- data[samples, , 
                          on=.(sample = sample_id)][batch_id != batch, 
                                                    median(beta.value, 
                                                           na.rm = TRUE),
                                                    by=.(feature)]$V1
    DF <- data.frame(median_batch - median_others, 
                     row.names = rownames(median_batch))
    colnames(DF) <- batch
    return(DF)
        # dif <- as.data.frame(matrix(NA, nrow = length(genes), ncol = 1))
        # rownames(dif) <- genes
        # for (i in seq_len(nrow(medianDif))) {
        #     j <- which(colnames(medianDif) == batch)
        #     ## barcode ids of all samples from one batch
        #     batchSamples <- samples$sample_id[samples$batch_id ==
        #                                           colnames(medianDif)[j]]
        #     ## barcode ids of all samples from all other batches
        #     otherSamples <- samples$sample_id[samples$batch_id !=
        #                                           colnames(medianDif)[j]]
        #     
        #     
        #     ## save medians for one gene in dif list
        #     dif[i, 1] <-
        #         abs(median(as.numeric(data[rownames(medianDif)[i],
        #                                    as.character(batchSamples)]), 
        #                    na.rm = TRUE) 
        #             - median(as.numeric(data[rownames(medianDif)[i],
        #                                        as.character(otherSamples)]), 
        #                        na.rm = TRUE))
        # }
        # return(dif)
    }