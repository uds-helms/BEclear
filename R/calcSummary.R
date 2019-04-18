#' calcSummary
#'
#' @aliases calcSummary
#'
#' @title Summarize median comparison and p-value calculation results
#'
#' @seealso \code{\link{calcBatchEffects}}
#' @seealso \code{\link{correctBatchEffect}}
#'
#' @description  Summarizes the results of the \code{\link{calcBatchEffects}} 
#' function
#'
#' @details All genes with a median comparison value >= 0.05 and a p-value
#' of <= 0.01 are summarized into a data.frame. These genes are assumed to
#' contain a batch effect
#'
#' @param medians a matrix containing median difference values calculated by
#' the \code{\link{calcBatchEffects}} function. For further details look at the
#' documentation of this function.
#' @param pvalues a matrix containing p-values calculated by the
#' \code{\link{calcBatchEffects}} function. For further details look at the
#' documentation of this function.
#' @param mediansTreshold the threshold above or equal median values are regarded 
#' as batch effected, when the criteria for p-values is also met.
#' @param pvaluesTreshold the threshold below or equal p-values are regarded as 
#' batch effected, when the criteria for medians is also met.
#'
#' @export calcSummary
#' @import futile.logger
#' @usage calcSummary(medians, pvalues, mediansTreshold, pvaluesTreshold)
#'
#' @return Null if there are no batch effects detected, else
#'  A data frame with the columns "gene" containing the gene name,
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
#' ex.data <- ex.data[31:90, 7:26]
#' ex.samples <- ex.samples[7:26, ]
#' 
#' ## Calculate the batch effects
#' batchEffects <- calcBatchEffects(data = ex.data, samples = ex.samples,
#' adjusted = TRUE, method = "fdr")
#' med <- batchEffects$med
#' pvals <- batchEffects$pval
#' 
#' ## Summarize p-values and median differences for batch affected genes
#' sum <- calcSummary(medians = med, pvalues = pvals)
calcSummary <- function(medians, pvalues, mediansTreshold = 0.05, 
                        pvaluesTreshold = 0.01) {
  ## build summary table of found genes
  flog.info("Generating a summary table")
  summaryTable <- as.data.frame(matrix(ncol = 4))
  colnames(summaryTable) <- c("gene", "batch", "median", "pvalue")
  counter <- 1
  for (i in seq_len(nrow(medians))) {
    for (j in seq_len(ncol(medians))) {
      if ((is.na(medians[i, j]) | (medians[i, j] >= mediansTreshold))
      & (pvalues[i, j] <= pvaluesTreshold)) {
        summaryTable[counter, "gene"] <- rownames(medians)[i]
        summaryTable[counter, "batch"] <- colnames(medians)[j]
        summaryTable[counter, "median"] <- medians[i, j]
        summaryTable[counter, "pvalue"] <- pvalues[i, j]
        counter <- counter + 1
      }
    }
  }
  if (counter == 1) {
    return(NULL)
  } else {
    return(summaryTable)
  }
}
