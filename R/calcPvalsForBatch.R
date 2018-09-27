#' calcPvalsForBatch
#'
#'@keywords internal
#'
#'@description  function to calculate p-values for every gene in a batch
#'
#'@importFrom stats ks.test
#'
#'@return the p-values for genes in a batch
calcPvalsForBatch <- function(batch, genes, pvalues, samples, data) {
    pvals <- as.data.frame(matrix(NA, nrow=length(genes), ncol=1))
    rownames(pvals) <- genes
    for(i in seq_len(nrow(pvalues))) {
        j <- which(colnames(pvalues) == batch)
        ## barcode ids of all samples from one batch
        batchSamples <- samples$sample_id[samples$batch_id == 
                                              colnames(pvalues)[j]]
        ## barcode ids of all samples from all other batches
        otherSamples <- samples$sample_id[samples$batch_id != 
                                              colnames(pvalues)[j]]
        
        x <- as.numeric(data[rownames(pvalues)[i], as.character(batchSamples)])
        if(all(is.na(x))){
            ## if all beta values of a gene in a batch are NA, the p-value is 
            ## set to zero
            pvals[i, 1] <- 0.0
        } else{
            y <- as.numeric(data[rownames(pvalues)[i], 
                                 as.character(otherSamples)])
            ## save pvalue for one gene in pvalue matrix
            ksResult <- ks.test(x, y)
            pvals[i, 1] <- ksResult$p.value
        }

    }
    return(pvals)
}