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
#' other batches and returns a p-value which defines if the distributions are 
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
#' @importFrom rlist list.cbind
#' @usage calcPvalues(data, samples, adjusted=TRUE, method="fdr", 
#' BPPARAM=bpparam())
#' 
#' @return a matrix containing p-values for all genes in all batches, the 
#' column names define the batch numbers, row names are the same gene names as 
#' contained in the input matrix. If there are only missing values present for
#' a gene in a batch, a p-values of 0 is returned
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

    flog.info(paste("Calculating the p-values for", samples[,unique(batch_id)]
                    , "batches"))
    
    result <- bplapply(samples[,unique(batch_id)], calcPvalsForBatch,
                       samples = samples, data = data, BPPARAM=BPPARAM)
    
    pvalues <- list.cbind(result)
    
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
