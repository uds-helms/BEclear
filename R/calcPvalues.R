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
#' @param data a \code{\link[data.table]{data.table}} with one column indicating
#' the sample, one the features and a value column indicating the beta value
#' or a matrix with rows as features and columns as samples
#' @param samples data frame with two columns, the first column has to contain 
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param adjusted should the p-values be adjusted or not, see "method" for 
#' available adjustment methods.
#' @param method adjustment method for p-value adjustment (if TRUE), default 
#' method is "false discovery rate adjustment", for other available methods 
#' see the description of the used standard R package \code{\link{p.adjust}}.
#' @param BPPARAM An instance of the 
#' \code{\link[BiocParallel]{BiocParallelParam-class}} that determines how to
#' parallelisation of the functions will be evaluated.
#' 
#' @export calcPvalues
#' @import futile.logger
#' @import data.table
#' @importFrom stats p.adjust
#' @importFrom methods is
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
#'  
#' 
#' pvals <- calcPvalues(data=ex.data, samples=ex.samples,method="fdr")

calcPvalues <- function(data, samples, adjusted=TRUE, method="fdr", 
                        BPPARAM=bpparam()) {
    
    if(!is(data, "data.table")){
        flog.info("Transforming matrix to data.table")
        data <- data.table(feature=as.character(rownames(data)), data)
        data <- melt(data = data, id.vars = "feature", variable.name = "sample", 
                   value.name = "beta.value", variable.factor = FALSE)
        setkey(data, "feature", "sample")
    }
    
    if(!is(samples, "data.table")){
        samples <- data.table(samples)
    }

    flog.info(paste("Calculating the p-values for", samples[,uniqueN(batch_id)] 
                    , "batches"))
    
    pvalues <- lapply(samples[,unique(batch_id)], calcPvalsForBatch,
                       samples = samples, data = data, BPPARAM=BPPARAM)
    
    flog.debug(paste("Binding", length(pvalues), "rows of p-values together"))
    pvalues <- do.call(cbind, pvalues)
    
    ## pvalue adjustment
    if (adjusted == TRUE) {
        flog.info("Adjusting p-values")
        pvalues <- t(apply(pvalues, 1, p.adjust, method=method))
    }
    
    return(pvalues)
}
