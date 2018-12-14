#' countValuesToPredict
#'
#' @aliases countValuesToPredict
#'
#' @seealso \code{\link{clearBEgenes}}
#' @seealso \code{\link{BEclear}}
#' @seealso \code{\link{correctBatchEffect}}
#'
#' @title Count NA entries in a matrix
#'
#' @description Simple function that counts all values in a matrix which are NA
#'
#' @return Returns a data frame with the number of NA entries per column. Since
#' the function is mainly written for the usage in batch effect correction of
#' DNA methylation data, the column names of the data frame are set to "sample"
#' and "num_pred_values". Nevertheless, the function can be used with any other
#' matrix containing anything but beta values.
#'
#' @param data any kind of matrix
#'
#' @export countValuesToPredict
#' @import futile.logger
#' @usage countValuesToPredict(data)
#'
#' @examples
#' ## Shortly running example. For a more realistic example that takes
#' ## some more time, run the same procedure with the full BEclearData
#' ## dataset
#' 
#' ## Whole procedure that has to be done to use this function.
#' data(BEclearData)
#' ex.data <- ex.data[31:90, 7:26]
#' ex.samples <- ex.samples[7:26, ]
#' 
#' # Calculate median difference values and p-values
#' library(data.table)
#' samples <- data.table(ex.samples)
#' data <- data.table(feature = rownames(ex.data), ex.data)
#' data <- melt(
#'   data = data, id.vars = "feature", variable.name = "sample",
#'   value.name = "beta.value"
#' )
#' setkey(data, "feature", "sample")
#' meds <- calcMedians(data = data, samples = samples)
#' pvals <- calcPvalues(data = data, samples = samples)
#' 
#' # Summarize p-values and median differences for batch affected genes
#' sum <- calcSummary(medians = meds, pvalues = pvals)
#' clearedMatrix <- clearBEgenes(data = ex.data, samples = ex.samples, summary = sum)
#' numberOfEntries <- countValuesToPredict(data = clearedMatrix)
countValuesToPredict <- function(data) {
  amount <- sum(is.na(data)) / (nrow(data) * ncol(data))
  flog.info(paste(sum(is.na(data)), "values (", amount * 100,
    "% of the data) set to NA",
    sep = " "
  ))
  numberPredValues <- matrix(ncol = 2, nrow = ncol(data))
  numberPredValues <- as.data.frame(numberPredValues)
  colnames(numberPredValues) <- c("sample", "num_pred_values")
  numberPredValues[, 1] <- colnames(data)
  numberPredValues[, 2] <- 0
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(ncol(data))) {
      if (is.na(data[i, j])) {
        numberPredValues[numberPredValues$sample == j, 2] <- numberPredValues[numberPredValues$sample == j, 2] + 1
      }
    }
  }
  return(numberPredValues = numberPredValues)
}
