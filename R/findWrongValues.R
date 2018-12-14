#' findWrongValues
#'
#' @aliases findWrongValues
#'
#' @seealso \code{\link{replaceWrongValues}}
#' @seealso \code{\link{correctBatchEffect}}
#'
#' @title Find DNA methylation values out of the boundaries
#'
#' @description A method which lists values below 0 or beyond 1 contained in the
#' input matrix. The wrong entries are stored in a data.frame together with the
#' corresponding row and column position of the matrix. Note that this method is
#' especially designed for DNA methylation data.
#'
#' @details Note that this method is especially designed to run after the batch
#' effect correction of DNA methylation data, e.g. with the
#' \code{\link{BEclear}} method. It can happen, that the predicted values are
#' lying slightly below the lower bound of 0 or beyond the upper bound of 1.
#' This method finds these inaccurately predicted entries. Another method called
#' \code{\link{replaceWrongValues}} replaces these values either by 0 or 1,
#' respectively.
#'
#' @param data any matrix filled with values that normally should be bounded
#' between 0 and 1.
#'
#' @export findWrongValues
#' @import futile.logger
#' @usage findWrongValues(data)
#'
#' @return A data frame containing the columns "level", "row", "col" and "value"
#' defining if the wrong value is below 0 or beyond 1 (level = 0 or level = 1),
#' the row position and the column position in the input matrix and the wrong
#' value itself, respectively.
#'
#' @examples
#' data(BEclearCorrected)
#' # Find wrongly predicted values
#' wrongEntries <- findWrongValues(data = ex.corrected.data)
findWrongValues <- function(data) {
  ## find entries > 1 or < 0 in the matrix, store wrong entries in a
  ## data.frame
  wrongEntries <- as.data.frame(matrix(nrow = 0, ncol = 4))
  colnames(wrongEntries) <- c("level", "row", "col", "value")
  counter <- 1
  for (i in seq_len(nrow(data))) {
    for (j in seq_len(ncol(data))) {
      if (data[i, j] > 1) {
        wrongEntries[counter, "level"] <- 1
        wrongEntries[counter, "row"] <- i
        wrongEntries[counter, "col"] <- j
        wrongEntries[counter, "value"] <- data[i, j]
        counter <- counter + 1
      }
      if (data[i, j] < 0) {
        wrongEntries[counter, "level"] <- 0
        wrongEntries[counter, "row"] <- i
        wrongEntries[counter, "col"] <- j
        wrongEntries[counter, "value"] <- data[i, j]
        counter <- counter + 1
      }
    }
  }
  flog.info(paste(counter - 1, " wrong values found", sep = ""))
  return(wrongEntries)
}
