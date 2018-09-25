#' replaceWrongValues
#' 
#' @aliases replaceWrongValues
#' 
#' @title Replace DNA methylation values out of the boundaries
#' 
#' @seealso \code{\link{findWrongValues}}
#' @seealso \code{\link{correctBatchEffect}}
#' 
#' @description A method which replaces values below 0 or beyond 1 contained in
#' the input matrix. These wrong entries are replaced by 0 or 1, respectively.
#' Note that this method is especially designed for DNA methylation data.
#' 
#' @details Note that this method is especially designed to run after the batch 
#' effect correction of DNA methylation data, e.g. with the 
#' \code{\link{imputeMissingData}} method. It can happen, that the predicted 
#' values are lying slightly below the lower bound of 0 or beyond the upper 
#' bound of 1. This method finds these inaccurately predicted entries. Another 
#' method called \code{\link{replaceWrongValues}} replaces these values either 
#' by 0 or 1, respectively. Another method called \code{\link{findWrongValues}} 
#' returns a list of existing wrong values and can be run before the 
#' replacement.
#' 
#' @param data any matrix filled with values that normally should be bounded 
#' between 0 and 1.
#' 
#' @export replaceWrongValues
#' @import futile.logger
#' @usage replaceWrongValues(data)
#' 
#' @return Returns the input matrix with every value previously below 0 changed 
#' to 0 and every value previously beyond 1 changed to 1.
#' 
#' @examples 
#' data(BEclearCorrected)
#' # Replace wrongly predicted values
#' corrected <- replaceWrongValues(data=ex.corrected.data)

replaceWrongValues <- function(data) {
  flog.info("Replacing values below 0 or above 1:")
    counter <- 0
    for(i in seq_len(nrow(data))) {
        for(j in seq_len(ncol(data))) {
            if(data[i, j] > 1) {
                data[i, j] <- 1
                counter <- counter + 1
            }
            if(data[i, j] < 0) {
                data[i, j] <- 0
                counter <- counter + 1
            }
        }
    }
    flog.info(paste(counter, "values replaced"))
    
    return(data)
}

