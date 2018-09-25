#' calcPvalues
#' 
#' @aliases calcPvalues
#' 
#' @title Calculate p-values
#' 
#' @seealso \code{\link{ks.test}}
#' @seealso \code{\link{p.adjust}}
#' @seealso \code{\link{correctBatchEffect}}
#' 
#' @description Compares the distribution of all beta values corresponding to 
#' one batch with the distribution of all beta values corresponding to all 
#' other batches and retuns a p-value which defines if the distributions are 
#' the same or not. Standard two sided Kolmogorov-Smirnov test is used to 
#' calculate the (adjusted) p-values.
#' 
#' 
#' @param data any matrix filled with beta values, column names have to be 
#' sample_ids corresponding to the ids listed in "samples", row names have to
#' be gene names.
#' @param samples data frame with two columns, the first column has to contain 
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param adjusted should the p-values be adjusted or not, see "method" for 
#' available adjustment methods.
#' @param method adjustment method for p-value adjustment (if TRUE), default 
#' method is "false discovery rate adjustment", for other available methods 
#' see the description of the used standard R package \code{\link{p.adjust}}.
#' @param BPPARAM An instance of the 
#' \code{\link[BiocParallel]{BiocParallelParam}} that determines how to
#' parallelisation of the functions will be evaluated.
#' 
#' @export calcPvalues
#' @import BiocParallel
#' @import futile.logger
#' @importFrom stats p.adjust
#' @usage calcPvalues(data, samples, adjusted=TRUE, method="fdr", 
#' BPPARAM=bpparam())
#' 
#' @return a matrix containing p-values for all genes in all batches, the 
#' column names define the batch numbers, row names are the same gene names as 
#' contained in the input matrix.
#' 
#' @examples 
#' ## Shortly running example. For a more realistic example that takes 
#' ## some more time, run the same procedure with the full BEclearData
#' ## dataset.
#' 
#' ## Calculate fdr-adjusted p-values in non-parallel mode
#' data(BEclearData)
#' ex.data <- ex.data[31:90,7:26]
#' ex.samples <- ex.samples[7:26,] 
#' pvals <- calcPvalues(data=ex.data, samples=ex.samples,method="fdr")

calcPvalues <- function(data, samples, adjusted=TRUE, method="fdr", 
                        BPPARAM=bpparam()) {
    ## get batch numbers
    batches <- unique(samples$batch_id)
    flog.info(paste("Calculating the p-values for", length(batches), "batches"))
    
    ## get genes
    genes <- rownames(data)
    
    ## construct data.frames filled with NA, one row per gene, one column per 
    ## batch
    pvalues <- as.data.frame(matrix(NA, nrow=length(genes),
                                    ncol=length(batches)))
    rownames(pvalues) <- genes
    colnames(pvalues) <- batches
    
    result <- bplapply(batches, calcPvalsForBatch, genes = genes, 
                       pvalues = pvalues, samples = samples, data = data,
                       BPPARAM=BPPARAM)
    
    ## fill pvalue matrix from result
    result <- unlist(result)
    counter <- 1
    for(j in seq_len(ncol(pvalues))) {
        for(i in seq_len(nrow(pvalues))) {
            pvalues[i, j] <- result[counter]
            counter <- counter + 1
        }
    }
    remove(result, counter, i, j)
    
    ## pvalue adjustment
    if (adjusted == TRUE) {
        pvaluesAdjusted <- pvalues
        ## p-value adjustment with false discovery rate
        for(i in seq_len(nrow(pvaluesAdjusted))) {
            pvaluesAdjusted[i, ] <- p.adjust(pvaluesAdjusted[i, ],
                                             method=method)
        }
        pvalues <- pvaluesAdjusted
    }
    
    return(pvalues)
}
