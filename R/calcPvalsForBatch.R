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
        ## save pvalue for one gene in pvalue matrix
        ksResult <- ks.test(as.numeric(data[rownames(pvalues)[i],
                                as.character(batchSamples)]),
                as.numeric(data[rownames(pvalues)[i],
                                as.character(otherSamples)]))
        pvals[i, 1] <- ksResult$p.value
    }
    return(pvals)
}