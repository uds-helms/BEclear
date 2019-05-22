#' clearBEgenes
#'
#' @aliases clearBEgenes
#'
#' @seealso \code{\link{calcBatchEffects}}
#' @seealso \code{\link{calcSummary}}
#' @seealso \code{\link{correctBatchEffect}}
#'
#' @title Prepare a data matrix for the BEclear function
#'
#' @description A function that simply sets all values to NA which were
#' previously found by median value comparison and p-value calculation and are
#' stored in a summary. The summary defines which values in the data matrix are
#' set to NA.
#'
#' @details All entries belonging to genes stated in the summary are set to NA
#' for the corresponding batches in the data matrix. Please look at the
#' descriptions of \code{\link{calcBatchEffects}}  for
#' more detailed information about the data which should be contained in the
#' summary data.frame.
#'
#' @param data any matrix filled with beta values, column names have to be
#' sample_ids corresponding to the ids listed in "samples", row names have to be
#' gene names.
#' @param samples data frame with two columns, the first column has to contain
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param summary a summary data.frame containing the columns "gene", "batch",
#' "median" and "p-value" and consists of all genes which were found in the
#' median and p-value calculations, see \code{\link{calcSummary}} function for
#' more details.
#'
#' @export clearBEgenes
#' @import futile.logger
#' @usage clearBEgenes(data, samples, summary)
#'
#' @return A data matrix with the same dimensions as well as the same column
#' and row names as the input data matrix is returned, all entries which are
#' defined in the summary are now set to NA.
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
#' meds <- batchEffects$med
#' pvals <- batchEffects$pval
#' 
#' ## Summarize p-values and median differences for batch affected genes
#' sum <- calcSummary(medians = meds, pvalues = pvals)
#' 
#' ## Set values for summarized BEgenes to NA
#' clearedMatrix <- clearBEgenes(data = ex.data, samples = ex.samples, summary = sum)
clearBEgenes <- function(data, samples, summary) {
  ## set beta values in data to NA for all found genes
  for (i in seq_len(nrow(summary)))
  {
    data[summary$gene[i], as.character(
      samples$sample_id[samples$batch_id == summary$batch_id[i]]
    )] <- NA
  }
  amount <- sum(is.na(data)) / (nrow(data) * ncol(data))
  flog.info("Removing values with batch effect:")
  flog.info(paste(sum(is.na(data)), "values (", amount * 100,
    "% of the data) set to NA",
    sep = " "
  ))
  return(data)
}
