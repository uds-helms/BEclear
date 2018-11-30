#' calcMedians
#' 
#' @aliases calcMedians
#' 
#' @title Calculate median difference values
#' 
#' @seealso \code{\link{correctBatchEffect}}
#' 
#' @description Compares the median value of all beta values belonging to one 
#' batch with the median value of all beta values belonging to all other 
#' batches. Returns a matrix containing this median difference value for every 
#' gene in every batch, columns define the batch numbers, rows the gene names.
#' 
#' 
#' @param data a \code{\link[data.table]{data.table}} with one column indicating
#' the sample, one the features and a value column indicating the beta value 
#' or a matrix with rows as features and columns as samples
#' @param samples data frame with two columns, the first column has to contain 
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param BPPARAM An instance of the 
#' \code{\link[BiocParallel]{BiocParallelParam-class}} that determines how to 
#' parallelisation of the functions will be evaluated.
#' 
#' @export calcMedians
#' @import futile.logger
#' @import data.table
#' @importFrom methods is
#' @usage calcMedians(data, samples, BPPARAM=bpparam())
#' 
#' @return a matrix containing median comparison values for all genes in all 
#' batches, the column names define the batch numbers, row names are the same 
#' gene names as contained in the input matrix.
#' 
#' @examples 
#' ## Shortly running example. For a more realistic example that takes
#' ## some more time, run the same procedure with the full BEclearData
#' ## dataset.
#' 
#' ## Calculate median comparison values in non-parallel mode
#' data(BEclearData)
#' ex.data <- ex.data[31:90,7:26]
#' ex.samples <- ex.samples[7:26,]
#' 
#' 
#' medians <- calcMedians(data=ex.data, samples=ex.samples)

calcMedians <- function(data, samples, BPPARAM=bpparam()) {
    .Deprecated("calcBatchEffects")
    
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

    flog.info(paste("Calculating the medians for", samples[,uniqueN(batch_id)],
                    "batches"))
    
    
    result <- lapply(samples[,unique(batch_id)], calcMediansForBatch, 
                        samples = samples, data = data, BPPARAM=BPPARAM)
    
    return(do.call(cbind, result))
}
