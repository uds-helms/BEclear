#' calcSummary
#' 
#' @aliases calcSummary
#' 
#' @title Summarize median comparison and p-value calculation results
#' 
#' @seealso \code{\link{calcMedians}}
#' @seealso \code{\link{calcPvalues}}
#' @seealso \code{\link{correctBatchEffect}}
#' 
#' @description Summarizes the results of the median comparison function 
#' \code{\link{calcMedians}} and the p-value calculation function 
#' \code{\link{calcPvalues}}. Should be used with the matrices originating from
#' these two functions.
#' 
#' @details All genes with a median comparison value >= 0.05 and a p-value 
#' of <= 0.01 are summarized into a data.frame. These genes are assumed to 
#' contain a batch effect
#' 
#' @param medians a matrix containing median difference values calculated by 
#' the \code{\link{calcMedians}} function. For further details look at the
#' documentation of this function.
#' @param pvalues a matrix containing p-values calculated by the 
#' \code{\link{calcPvalues}} function. For further details look at the 
#' documentation of this function.
#'
#' @export calcSummary
#' @import futile.logger
#' @usage calcSummary(medians, pvalues)
#' 
#' @return A data frame with the columns "gene" containing the gene name, 
#' "batch" containing the batch number from which the gene was found, "median" 
#' and "p-value" containing the calculated median difference values and the 
#' p-values, respectively.
#' 
#' @examples 
#' ## Shortly running example. For a more realistic example that takes
#' ## some more time, run the same procedure with the full BEclearData
#' ## dataset.
#' 
#' ## Whole procedure that has to be done to use this function.
#' data(BEclearData)
#' ex.data <- ex.data[31:90,7:26]
#' ex.samples <- ex.samples[7:26,]
#' 
#' # Calculates median difference values and p-values from the example data
#'library(data.table)
#'samples <- data.table(ex.samples)
#'data <- data.table(feature=rownames(ex.data), ex.data)
#'data <- melt(data = data, id.vars = "feature", variable.name = "sample", 
#'             value.name = "beta.value")
#'setkey(data, "feature", "sample")
#'med <- calcMedians(data=data, samples=samples)
#' pvals <- calcPvalues(data=data, samples=samples, adjusted=TRUE, 
#' method="fdr")
#' 
#' # Summarize p-values and median differences for batch affected genes 
#' sum <- calcSummary(medians=med, pvalues=pvals)

calcSummary <- function(medians, pvalues) {
    ## build summary table of found genes
    flog.info("Generating a summary table")
    summaryTable <- as.data.frame(matrix(ncol=4))
    colnames(summaryTable) <- c("gene", "batch", "median", "pvalue")
    counter <-1
    for(i in seq_len(nrow(medians))) {
        for(j in seq_len(ncol(medians))) {
            if((is.na(medians[i, j]) | (medians[i, j] >= 0.05)) 
                 & (pvalues[i, j] <= 0.01)) {
                summaryTable[counter, "gene"] <- rownames(medians)[i]
                summaryTable[counter, "batch"] <- colnames(medians)[j]
                summaryTable[counter, "median"] <- medians[i, j]
                summaryTable[counter, "pvalue"] <- pvalues[i, j]
                counter <- counter + 1
            }
        }
    }
    return(summaryTable)
}
