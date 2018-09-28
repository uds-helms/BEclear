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
#' @param data any matrix filled with beta values, column names have to be 
#' sample_ids corresponding to the ids listed in "samples", row names have to be
#' gene names.
#' @param samples data frame with two columns, the first column has to contain 
#' the sample numbers, the second column has to contain the corresponding batch
#' number. Colnames have to be named as "sample_id" and "batch_id".
#' @param BPPARAM An instance of the 
#' \code{\link[BiocParallel]{BiocParallelParam}} that determines how to 
#' parallelisation of the functions will be evaluated.
#' 
#' @export calcMedians
#' @import BiocParallel
#' @import futile.logger
#' @import data.table
#' @importFrom rlist list.cbind
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
#' library(data.table)
#' samples <- data.table(ex.samples)
#' data <- data.table(feature=rownames(ex.data), ex.data)
#' data <- melt(data = data, id.vars = "feature", variable.name = "sample", 
#' value.name = "beta.value")
#' setkey(data, "feature", "sample")
#' 
#' medians <- calcMedians(data=data, samples=samples)

calcMedians <- function(data, samples, BPPARAM=bpparam()) {

    flog.info(paste("Calculating the medians for", samples[,uniqueN(batch_id)],
                    "batches"))
    
    
    result <- bplapply(samples[,unique(batch_id)], calcMediansForBatch, 
                        samples = samples, data = data, BPPARAM=BPPARAM)
    
    # ## fill median matrix from result
    # result <- unlist(result)
    # counter <- 1
    # for(j in seq_len(ncol(medianDif))) {
    #     for(i in seq_len(nrow(medianDif))) {
    #         medianDif[i, j] <- result[counter]
    #         counter <- counter + 1
    #     }
    # }
    # remove(result, counter, i, j)
    
    return(list.cbind(result))
}
