#' findOutsideValues
#'
#' @aliases findOutsideValues
#'
#' @seealso \code{\link{replaceOutsideValues}}
#' @seealso \code{\link{correctBatchEffect}}
#'
#' @title Find DNA methylation values out of the boundaries
#'
#' @description A method which lists values below 0 or beyond 1 contained in the
#' input matrix. These entries are stored in a data.frame together with the
#' corresponding row and column position of the matrix. Note that this method is
#' especially designed for DNA methylation data.
#'
#' @details Note that this method is especially designed to run after the batch
#' effect correction of DNA methylation data, e.g. with the
#' \code{\link{BEclear}} method. It can happen, that the predicted values are
#' lying slightly below the lower bound of 0 or beyond the upper bound of 1.
#' This method finds these inaccurately predicted entries. Another method called
#' \code{\link{replaceOutsideValues}} replaces these values either by 0 or 1,
#' respectively.
#'
#' @param data any matrix filled with values that normally should be bounded
#' between 0 and 1.
#'
#' @export findOutsideValues
#' @import futile.logger
#' @usage findOutsideValues(data)
#'
#' @return A data frame containing the columns "level", "row", "col" and "value"
#' defining if the value is below 0 or beyond 1 (level = 0 or level = 1),
#' the row position and the column position in the input matrix and the
#' value itself, respectively.
#'
#' @examples
#' data(BEclearCorrected)
#' # Find predicted values outside of the boundaries
#' outsideEntries <- findOutsideValues(data = ex.corrected.data)
findOutsideValues <- function(data) {
  ## find entries > 1 or < 0 in the matrix, store these entries in a
  ## data.frame
  outsideEntries <- as.data.frame(matrix(nrow = 0, ncol = 4))
  colnames(outsideEntries) <- c("level", "row", "col", "value")
  counter <- 1
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(ncol(data))) {
      if (data[i, j] > 1) {
        outsideEntries[counter, "level"] <- 1
        outsideEntries[counter, "row"] <- i
        outsideEntries[counter, "col"] <- j
        outsideEntries[counter, "value"] <- data[i, j]
        counter <- counter + 1
      }
      if (data[i, j] < 0) {
        outsideEntries[counter, "level"] <- 0
        outsideEntries[counter, "row"] <- i
        outsideEntries[counter, "col"] <- j
        outsideEntries[counter, "value"] <- data[i, j]
        counter <- counter + 1
      }
    }
  }
  flog.info(paste(counter - 1, " values outside of the boundaries found", sep = ""))
  return(outsideEntries)
}
